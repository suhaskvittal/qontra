/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#include "graph/decoding_graph.h"

namespace qontra {
namespace graph {

#define N_COORD 3

typedef std::function<void(fp_t, std::vector<uint>, std::set<uint>)>
    error_callback_t;
typedef std::function<void(uint, std::array<fp_t, N_COORD>)>
    detector_callback_t;

void
read_detector_error_model(const stim::DetectorErrorModel&, 
        uint n_iter, uint& det_offset, 
        std::array<fp_t, N_COORD>& coord_offset,
        error_callback_t, detector_callback_t);

using namespace decoding;

void
DecodingGraph::build_distance_matrix() {
    dijkstra::ewf_t<vertex_t> w = [&] (vertex_t* v1, vertex_t* v2)
    {
        auto e = this->get_edge(v1, v2);
        if (e == nullptr) {
            return 1000000000.0;
        } else {
            return e->edge_weight;
        }
    };
    dijkstra::callback_t<vertex_t, matrix_entry_t> cb = 
    [&] (vertex_t* src,
        vertex_t* dst,
        const std::map<vertex_t*, fp_t>& dist,
        const std::map<vertex_t*, vertex_t*>& pred)
    {
        uint32_t length = 0;
        std::set<uint> frames;

        fp_t weight = dist.at(dst);
        if (weight < 1000) {
            auto curr = dst;
            while (curr != src) {
                auto next = pred.at(curr);
                if (curr == next) {
                    weight = 1000000000;
                    goto failed;
                }
                auto e = this->get_edge(next, curr);
                std::set<uint> new_frames;
                std::set_symmetric_difference(
                        frames.begin(), frames.end(),
                        e->frames.begin(), e->frames.end(),
                        std::inserter(new_frames, new_frames.end()));
                // Update data
                frames = new_frames;
                length++;
                curr = next;
            }
        }
failed:
        fp_t prob = pow(10, -weight);

        return (matrix_entry_t) {length, prob, weight, frames};
    };

    distance_matrix = create_distance_matrix(this, w, cb);
}

void
DecodingGraph::build_error_polynomial() {
    auto e0 = edges[0];
    poly_t pX{1 - e0->error_probability, e0->error_probability};
    // Multiply the polynomials together
    fp_t expectation = 0.0;
    for (uint i = 1; i < edges.size(); i++) {
        auto e = edges[i];
        poly_t a(pX.size()+1);  // p(X) * (1-e)
        poly_t b(pX.size()+1);  // p(X) * eX
        
        b[0] = 0;
        for (uint j = 0; j < pX.size(); j++) {
            a[j] += pX[j] * (1 - e->error_probability);
            b[j+1] += pX[j] * e->error_probability;
        }
        pX.clear(); pX.push_back(0);    // Clear and increment the size of pX
        for (uint j = 0; j < pX.size(); j++) {
            pX[j] = a[0] + b[0];
            if (i == edges.size()-1) {
                expectation += j * pX[j];
            }
        }
    }
    error_polynomial = pX;
    expected_errors = expectation;
}

bool 
DecodingGraph::update_state() {
    if (!__DecodingGraphParent::update_state()) return false;
    build_distance_matrix();
    build_error_polynomial();
    return true;
}

DecodingGraph
to_decoding_graph(const stim::Circuit& qec_circ) {
    DecodingGraph graph;

    stim::DetectorErrorModel dem = 
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
            qec_circ,
            true,  // decompose_errors
            true,  // fold loops
            false, // allow gauge detectors
            1.0,   // approx disjoint errors threshold
            false, // ignore decomposition failures
            false
        );
    // Create callbacks.
    detector_callback_t det_f =
        [&graph](uint det, std::array<fp_t, N_COORD> coords) 
        {
            vertex_t* v;
            if ((v=graph.get_vertex(det)) == nullptr) {
                v = new vertex_t;
                v->id = det;
                graph.add_vertex(v);
            }
            v->coords = coords;
        };
    error_callback_t err_f = 
        [&graph](fp_t prob, std::vector<uint> dets, std::set<uint> frames)
        {
            if (prob == 0 || dets.size() == 0 || dets.size() > 2) {
                return;  // Zero error probability -- not an edge.
            }
            
            if (dets.size() == 1) {
                // We are connecting to the boundary here.
                dets.push_back(BOUNDARY_INDEX);
            }
            // Now, we should only have two entries in det.
            vertex_t* v1;
            vertex_t* v2;
            if ((v1=graph.get_vertex(dets[0])) == nullptr) {
                vertex_t* v = new vertex_t;
                v->id = dets[0];
                graph.add_vertex(v);
                v1 = v;
            }
            if ((v2=graph.get_vertex(dets[1])) == nullptr) {
                vertex_t* v = new vertex_t;
                v->id = dets[1];
                graph.add_vertex(v);
                v2 = v;
            }
            auto e = graph.get_edge(v1, v2);
            if (e != nullptr) {
                fp_t old_prob = e->error_probability;
                std::set<uint> old_frames = e->frames;
                if (frames == old_frames) {
                    prob = prob * (1-old_prob) + old_prob * (1-prob);
                }
            } else {
                // Create new edge if it does not exist.
                e = new edge_t;
                e->src = (void*)v1;
                e->dst = (void*)v2;
                graph.add_edge(e);
            }
            fp_t edge_weight = (fp_t)log10((1-prob)/prob);
            e->edge_weight = edge_weight;
            e->error_probability = prob;
            e->frames = frames;
        };
    // Declare coord offset array.
    uint det_offset = 0;
    std::array<fp_t, N_COORD> coord_offset;
    coord_offset.fill(0);  // Zero initialize.
    // Use callbacks to build graph.
    read_detector_error_model(dem, 1, det_offset, coord_offset,
                                err_f, det_f);
    return graph;
}

void 
read_detector_error_model(
        const stim::DetectorErrorModel& dem, uint n_iter,
        uint& det_offset, std::array<fp_t, N_COORD>& coord_offset,
        error_callback_t err_f, detector_callback_t det_f) 
{
    while (n_iter--) {  // Need this to handle repeats.
        for (stim::DemInstruction inst : dem.instructions) {
            stim::DemInstructionType type = inst.type;
            if (type == stim::DemInstructionType::DEM_REPEAT_BLOCK) {
                // The targets for this instruction are
                // (1) number of repeats, and
                // (2) block number.
                uint n_repeats = (uint)inst.target_data[0].data;
                stim::DetectorErrorModel subblock = 
                    dem.blocks[inst.target_data[1].data];
                read_detector_error_model(subblock, n_repeats,
                           det_offset, coord_offset, err_f, det_f);
            } else if (type == stim::DemInstructionType::DEM_ERROR) {
                std::vector<uint> detectors;
                std::set<uint> frames;
                 
                fp_t e_prob = (fp_t)inst.arg_data[0];
                for (stim::DemTarget target : inst.target_data) {
                    if (target.is_relative_detector_id()) {
                        // This is a detector, add it to the list.
                        detectors.push_back(
                                (uint)target.data + det_offset);
                    } else if (target.is_observable_id()) {
                        frames.insert(target.data); 
                    } else if (target.is_separator()) {
                        // This is just due to decomposition.
                        // Handle each part of the decomposition
                        // separately.
                        err_f(e_prob, detectors, frames);
                        // Clear detectors and frames.
                        // We have already done the callback.
                        detectors.clear();
                        frames.clear();
                    }
                }
                // Handle last error.
                err_f(e_prob, detectors, frames);
            } else if (type == stim::DemInstructionType::DEM_SHIFT_DETECTORS) {
                det_offset += inst.target_data[0].data;
                uint k = 0;
                for (double a : inst.arg_data) {
                    coord_offset[k++] += (fp_t)a;
                }
            } else if (type == stim::DemInstructionType::DEM_DETECTOR) {
                // Compute coordinates.
                std::array<fp_t, N_COORD> coords(coord_offset);
                uint k = 0;
                for (double a : inst.arg_data) {
                    coords[k++] += (fp_t)a; 
                }
                // Now go through all declared detectors.
                for (stim::DemTarget target : inst.target_data) {
                    det_f(target.data + det_offset, coords);
                }
            }
        }
    }
}

}   // graph
}   // qontra
