/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#include "qontra/graph/decoding_graph.h"

#include <vtils/set_algebra.h>

namespace qontra {
namespace graph {

typedef std::function<void(fp_t, std::vector<uint64_t>, std::set<uint64_t>)>
    error_callback_t;
typedef std::function<void(uint64_t, std::array<fp_t, N_COORD>)>
    detector_callback_t;

void
read_detector_error_model(
        const stim::DetectorErrorModel&, 
        uint64_t n_iter,
        uint64_t& det_offset, 
        std::array<fp_t, N_COORD>& coord_offset,
        error_callback_t,
        detector_callback_t);

using namespace decoding;

//
// DecodingGraph Definitions
//

DecodingGraph::DecodingGraph(Mode m)
    :Graph(),
    distance_matrix(),
    error_polynomial(),
    expected_errors(),
    mode(m)
{
    std::array<fp_t, N_COORD> boundary_coords;
    boundary_coords.fill(-1);

    sptr<decoding::vertex_t> boundary = std::make_shared<decoding::vertex_t>();
    boundary->id = BOUNDARY_INDEX;
    boundary->coords = boundary_coords;
    add_vertex(boundary);
}

DecodingGraph::DecodingGraph(const DecodingGraph& other)
    :Graph(other),
    distance_matrix(other.distance_matrix),
    error_polynomial(other.error_polynomial),
    expected_errors(other.expected_errors),
    mode(other.mode)
{}

DecodingGraph&
DecodingGraph::operator=(const DecodingGraph& other) {
    Graph::operator=(other);
    distance_matrix = other.distance_matrix;
    error_polynomial = other.error_polynomial;
    expected_errors = other.expected_errors;
    mode = other.mode;
    return *this;
}

DecodingGraph::matrix_entry_t
DecodingGraph::dijkstra_cb(sptr<vertex_t> src,
                            sptr<vertex_t> dst,
                            const std::map<sptr<vertex_t>, fp_t>& dist,
                            const std::map<sptr<vertex_t>, sptr<vertex_t>>& pred)
{
    uint32_t length = 0;
    std::set<uint64_t> frames;

    fp_t weight = dist.at(dst);
    bool found_boundary = false;
    std::vector<sptr<vertex_t>> path;
    if (weight < 1000) {
        auto curr = dst;
        while (curr != src) {
            auto next = pred.at(curr);
            if (curr == next) {
                weight = 1000000000;
                path.clear();
                found_boundary = false;
                goto failed;
            }
            auto e = get_edge(next, curr);
            // Update data
            frames ^= e->frames;
            length++;
            path.push_back(curr);
            found_boundary |= (is_boundary(curr) || is_colored_boundary(curr)) && (curr != dst);
            curr = next;
        }
        path.push_back(src);
    }
failed:
    fp_t prob = pow(10, -weight);
    return (matrix_entry_t) {length, prob, weight, frames, path, found_boundary};
}

void
DecodingGraph::build_error_polynomial() {
    sptr<edge_t> e0 = edges[0];
    poly_t pX{1 - e0->error_probability, e0->error_probability};
    // Multiply the polynomials together
    fp_t expectation = 0.0;
    for (size_t i = 1; i < edges.size(); i++) {
        sptr<edge_t> e = edges[i];
        poly_t a(pX.size()+1);  // p(X) * (1-e)
        poly_t b(pX.size()+1);  // p(X) * eX
        
        b[0] = 0;
        for (size_t j = 0; j < pX.size(); j++) {
            a[j] += pX[j] * (1 - e->error_probability);
            b[j+1] += pX[j] * e->error_probability;
        }
        pX.clear(); pX.push_back(0);    // Clear and increment the size of pX
        for (size_t j = 0; j < pX.size(); j++) {
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
    if (mode != Mode::LOW_MEMORY) {
        build_distance_matrix();
    }
    build_error_polynomial();
    return true;
}

DecodingGraph
to_decoding_graph(const stim::Circuit& circuit, DecodingGraph::Mode mode) {
    if (mode == DecodingGraph::Mode::DO_NOT_BUILD)  return DecodingGraph();
    DecodingGraph graph(mode);

    stim::DetectorErrorModel dem = 
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
            circuit,
            true,  // decompose_errors
            true,  // fold loops
            false, // allow gauge detectors
            1.0,   // approx disjoint errors threshold
            false, // ignore decomposition failures
            false
        );
    // Create callbacks.
    detector_callback_t det_f =
        [&graph](uint64_t det, std::array<fp_t, N_COORD> coords) 
        {
            sptr<vertex_t> v = graph.get_vertex(det);
            if (v == nullptr) {
                v = std::make_shared<vertex_t>();
                v->id = det;
                graph.add_vertex(v);
            }
            v->coords = coords;
        };
    error_callback_t err_f = 
        [&](fp_t prob, std::vector<uint64_t> dets, std::set<uint64_t> frames)
        {
            if (prob == 0 || dets.size() == 0 || dets.size() > 2) {
                return;  // Zero error probability -- not an edge.
            }
            
            if (dets.size() == 1) {
                // We are connecting to the boundary here.
                dets.push_back(BOUNDARY_INDEX);
            }
            // Now, we should only have two entries in det.
            sptr<vertex_t> v1 = graph.get_vertex(dets[0]);
            sptr<vertex_t> v2 = graph.get_vertex(dets[1]);
            if (v1 == nullptr) {
                sptr<vertex_t> v = std::make_shared<vertex_t>();
                v->id = dets[0];
                graph.add_vertex(v);
                v1 = v;
            }
            if (v2 == nullptr) {
                sptr<vertex_t> v = std::make_shared<vertex_t>();
                v->id = dets[1];
                graph.add_vertex(v);
                v2 = v;
            }
            auto e = graph.get_edge(v1, v2);
            if (e != nullptr) {
                fp_t old_prob = e->error_probability;
                std::set<uint64_t> old_frames = e->frames;
                if (frames == old_frames) {
                    prob = prob * (1-old_prob) + old_prob * (1-prob);
                }
            } else {
                // Create new edge if it does not exist.
                e = std::make_shared<edge_t>();
                e->src = std::static_pointer_cast<void>(v1);
                e->dst = std::static_pointer_cast<void>(v2);
                graph.add_edge(e);
            }
            fp_t edge_weight = (fp_t)log10((1-prob)/prob);
            e->edge_weight = edge_weight;
            e->error_probability = prob;
            e->frames = frames;
        };
    // Declare coord offset array.
    uint64_t det_offset = 0;
    std::array<fp_t, N_COORD> coord_offset;
    coord_offset.fill(0);  // Zero initialize.
    // Use callbacks to build graph.
    read_detector_error_model(dem, 1, det_offset, coord_offset, err_f, det_f);
    return graph;
}

void 
read_detector_error_model(
        const stim::DetectorErrorModel& dem,
        uint64_t n_iter,
        uint64_t& det_offset,
        std::array<fp_t, N_COORD>& coord_offset,
        error_callback_t err_f,
        detector_callback_t det_f) 
{
    while (n_iter--) {  // Need this to handle repeats.
        for (stim::DemInstruction inst : dem.instructions) {
            stim::DemInstructionType type = inst.type;
            if (type == stim::DemInstructionType::DEM_REPEAT_BLOCK) {
                // The targets for this instruction are
                // (1) number of repeats, and
                // (2) block number.
                size_t n_repeats = static_cast<size_t>(inst.target_data[0].data);
                stim::DetectorErrorModel subblock = dem.blocks[inst.target_data[1].data];
                read_detector_error_model(
                        subblock, n_repeats, det_offset, coord_offset, err_f, det_f);
            } else if (type == stim::DemInstructionType::DEM_ERROR) {
                std::vector<uint64_t> detectors;
                std::set<uint64_t> frames;
                 
                fp_t e_prob = static_cast<fp_t>(inst.arg_data[0]);
                for (stim::DemTarget target : inst.target_data) {
                    if (target.is_relative_detector_id()) {
                        // This is a detector, add it to the list.
                        detectors.push_back(
                                static_cast<uint64_t>(target.data + det_offset));
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
                size_t k = 0;
                for (double a : inst.arg_data) {
                    coord_offset[k++] += (fp_t)a;
                }
            } else if (type == stim::DemInstructionType::DEM_DETECTOR) {
                // Compute coordinates.
                std::array<fp_t, N_COORD> coords(coord_offset);
                size_t k = 0;
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
