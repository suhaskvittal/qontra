/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#include "graph/decoding_graph.h"

namespace qontra {
namespace graph {

typedef std::function<void(fp_t, std::vector<uint64_t>, std::set<uint>)>
    error_callback_t;
typedef std::function<void(uint, std::array<fp_t, N_COORD>)>
    detector_callback_t;

void
read_detector_error_model(
        const stim::DetectorErrorModel&, 
        uint n_iter,
        uint& det_offset, 
        std::array<fp_t, N_COORD>& coord_offset,
        error_callback_t,
        detector_callback_t);

using namespace decoding;

std::string
print_v(sptr<colored_vertex_t> v) {
    std::string s;
    if (is_colored_boundary(v)) s += "B";
    else                s += std::to_string(v->id);
    s += "[" + v->color + "]";
    return s;
}

//
// DecodingGraph Definitions
//

void
DecodingGraph::setup_flagged_decoding_graph(
        const std::vector<sptr<vertex_t>>& detectors,
        const std::vector<flag_edge_t>& flag_edges)
{
    const fp_t EPS = 1e-8;

    build_flagged_decoding_graph();
    // Update edge weights in the graph.
    for (flag_edge_t fe : flag_edges) {
        auto src = std::get<0>(fe),
             thru = std::get<1>(fe),
             dst = std::get<2>(fe);
        auto e1 = flagged_decoding_graph->get_edge(src, thru);
        auto e2 = flagged_decoding_graph->get_edge(thru, dst);
        if (e1 == nullptr || e2 == nullptr) continue;
        e1->edge_weight = 1.0;
        e2->edge_weight = 1.0;
    }

    // Now, redo dijkstras.
    ewf_t<vertex_t> wf = [&] (sptr<vertex_t> v1, sptr<vertex_t> v2) {
        return flagged_decoding_graph->_ewf(v1, v2);
    };
    for (uint i = 0; i < detectors.size(); i++) {
        auto v = detectors[i];
        std::map<sptr<vertex_t>, fp_t> dist;
        std::map<sptr<vertex_t>, sptr<vertex_t>> pred;
        distance::dijkstra(flagged_decoding_graph.get(), v, dist, pred, wf);
        // Update distance matrix.
        for (uint j = 0; j < detectors.size(); j++) {
            if (i == j) continue;
            auto w = detectors[j];
            auto new_entry = flagged_decoding_graph->_dijkstra_cb(v, w, dist, pred);
            tlm::put(flagged_decoding_graph->distance_matrix, v, w, new_entry);
        }
    }
    graph_has_changed = false;  // Avoid later updates.
}

fp_t
DecodingGraph::_ewf(sptr<vertex_t> v1, sptr<vertex_t> v2) {
    auto e = get_edge(v1, v2);
    return e->edge_weight;
}

DecodingGraph::matrix_entry_t
DecodingGraph::_dijkstra_cb(sptr<vertex_t> src,
                            sptr<vertex_t> dst,
                            const std::map<sptr<vertex_t>, fp_t>& dist,
                            const std::map<sptr<vertex_t>, sptr<vertex_t>>& pred)
{
    uint32_t length = 0;
    std::set<uint> frames;

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
            std::set<uint> new_frames;
            std::set_symmetric_difference(
                    frames.begin(), frames.end(),
                    e->frames.begin(), e->frames.end(),
                    std::inserter(new_frames, new_frames.end()));
            // Update data
            frames = new_frames;
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
DecodingGraph::build_distance_matrix() {
    ewf_t<vertex_t> wf = [&] (sptr<vertex_t> v1, sptr<vertex_t> v2) {
        return this->_ewf(v1, v2);
    };
    distance::callback_t<vertex_t, matrix_entry_t> d_cb = 
    [&] (sptr<vertex_t> v1, 
            sptr<vertex_t> v2,
            const std::map<sptr<vertex_t>, fp_t>& dist,
            const std::map<sptr<vertex_t>, sptr<vertex_t>>& pred)
    {
        return this->_dijkstra_cb(v1, v2, dist, pred);
    };
    distance_matrix = distance::create_distance_matrix(this, wf, d_cb);
}

void
DecodingGraph::build_flagged_decoding_graph() {
    flagged_decoding_graph = std::make_unique<DecodingGraph>();
    // Copy vertices over, but make new edge pointers.
    for (auto v : get_vertices()) flagged_decoding_graph->add_vertex(v);
    for (auto e : get_edges()) {
        sptr<edge_t> _e = std::make_shared<edge_t>();
        _e->src = e->src;
        _e->dst = e->dst;
        _e->edge_weight = e->edge_weight;
        _e->error_probability = e->error_probability;
        _e->frames = e->frames;
        flagged_decoding_graph->add_edge(_e);
    }
    // Copy distance matrix over.
    flagged_decoding_graph->distance_matrix = distance_matrix;
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
        [&](fp_t prob, std::vector<uint64_t> dets, std::set<uint> frames)
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
                std::set<uint> old_frames = e->frames;
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
    uint det_offset = 0;
    std::array<fp_t, N_COORD> coord_offset;
    coord_offset.fill(0);  // Zero initialize.
    // Use callbacks to build graph.
    read_detector_error_model(dem, 1, det_offset, coord_offset, err_f, det_f);
    return graph;
}

//
// ColoredDecodingGraph Definitions
//

face_t
make_face(sptr<colored_vertex_t> v1, sptr<colored_vertex_t> v2, sptr<colored_vertex_t> v3) {
    std::vector<sptr<colored_vertex_t>> vertices{v1, v2, v3};
    std::sort(vertices.begin(), vertices.end());
    return std::make_tuple(vertices[0], vertices[1], vertices[2]);
}

ColoredDecodingGraph::ColoredDecodingGraph(DecodingGraph::Mode mode)
    :Graph(),
    restricted_color_map(),
    restricted_graphs(),
    face_frame_map()
{
    restricted_color_map["rg"] = 0;
    restricted_color_map["gr"] = 0;
    restricted_color_map["rb"] = 1;
    restricted_color_map["br"] = 1;
    restricted_color_map["gb"] = 2;
    restricted_color_map["bg"] = 2;

    // Modify each decoding graph by introducing two boundaries -- one for each color.
    // Delete their existing boundaries.
    for (uint i = 0; i < 3; i++) {
        DecodingGraph gr(mode);
        auto boundary = gr.get_vertex(BOUNDARY_INDEX);
        gr.delete_vertex(boundary);
        restricted_graphs[i] = gr;
    }
    sptr<colored_vertex_t> bred = std::make_shared<colored_vertex_t>();
    sptr<colored_vertex_t> bgreen = std::make_shared<colored_vertex_t>();
    sptr<colored_vertex_t> bblue = std::make_shared<colored_vertex_t>();

    bred->id = RED_BOUNDARY_INDEX;
    bred->color = "r";
    bgreen->id = GREEN_BOUNDARY_INDEX;
    bgreen->color = "g";
    bblue->id = BLUE_BOUNDARY_INDEX;
    bblue->color = "b";

    add_vertex(bred);
    add_vertex(bgreen);
    add_vertex(bblue);
    // Add edges between each boundary.
    sptr<colored_edge_t> erg = std::make_shared<colored_edge_t>();
    sptr<colored_edge_t> erb = std::make_shared<colored_edge_t>();
    sptr<colored_edge_t> egb = std::make_shared<colored_edge_t>();

    erg->src = std::static_pointer_cast<void>(bred);
    erg->dst = std::static_pointer_cast<void>(bgreen);
    erg->is_undirected = true;
    erg->edge_weight = 1e-3;
    erg->error_probability = 1.0;

    erb->src = std::static_pointer_cast<void>(bred);
    erb->dst = std::static_pointer_cast<void>(bblue);
    erb->is_undirected = true;
    erb->edge_weight = 1e-3;
    erb->error_probability = 1.0;

    egb->src = std::static_pointer_cast<void>(bgreen);
    egb->dst = std::static_pointer_cast<void>(bblue);
    egb->is_undirected = true;
    egb->edge_weight = 1e-3;
    egb->error_probability = 1.0;

    add_edge(erg);
    add_edge(erb);
    add_edge(egb);
}

bool
ColoredDecodingGraph::add_vertex(sptr<colored_vertex_t> v) {
    if (!__ColoredDecodingGraphParent::add_vertex(v))   return false;
    for (auto& pair : restricted_color_map) {
        std::string r = pair.first;
        if (r[0] != v->color[0])    continue;
        DecodingGraph& gr = restricted_graphs[pair.second];
        if (!gr.add_vertex(v)) {
            __ColoredDecodingGraphParent::delete_vertex(v);
            return false;
        }
    }
    return true;
}

bool
ColoredDecodingGraph::add_edge(sptr<colored_edge_t> e) {
    if (!__ColoredDecodingGraphParent::add_edge(e)) return false;
    sptr<colored_vertex_t> src = std::reinterpret_pointer_cast<colored_vertex_t>(e->src);
    sptr<colored_vertex_t> dst = std::reinterpret_pointer_cast<colored_vertex_t>(e->dst);

    if (src->color == dst->color) {
        for (auto& pair : restricted_color_map) {
            std::string r = pair.first;
            DecodingGraph& gr = restricted_graphs[pair.second];
            if (r[0] == src->color[0] && !gr.add_edge(e)) {
                __ColoredDecodingGraphParent::delete_edge(e);
                return false;
            }
        }
    } else {
        std::string lat = src->color + dst->color;
        DecodingGraph& gr = restricted_graphs[restricted_color_map[lat]];
        if (!gr.add_edge(e)) {
            __ColoredDecodingGraphParent::delete_edge(e);
            return false;
        }
    }
    return true;
}

void
ColoredDecodingGraph::delete_vertex(sptr<colored_vertex_t> v) {
    // Delete the vertex in the corresponding decoding graphs.
    for (auto pair : restricted_color_map) {
        std::string r = pair.first;
        if (r[0] != v->color[0])    continue;
        DecodingGraph& gr = restricted_graphs[pair.second];
        gr.delete_vertex(v);
    }
    // Finally, delete the final reference to the vertex.
    __ColoredDecodingGraphParent::delete_vertex(v);
}

void
ColoredDecodingGraph::delete_edge(sptr<colored_edge_t> e) {
    sptr<colored_vertex_t> src = std::reinterpret_pointer_cast<colored_vertex_t>(e->src);
    sptr<colored_vertex_t> dst = std::reinterpret_pointer_cast<colored_vertex_t>(e->dst);

    if (src->color == dst->color) {
        for (auto& pair : restricted_color_map) {
            std::string r = pair.first;
            DecodingGraph& gr = restricted_graphs[pair.second];
            if (r[0] == src->color[0]) {
                gr.delete_edge(e);
            }
        }
    } else {
        std::string lat = src->color + dst->color;
        DecodingGraph& gr = restricted_graphs[restricted_color_map[lat]];
        gr.delete_edge(e);
    }
    __ColoredDecodingGraphParent::delete_edge(e);
}

bool
ColoredDecodingGraph::are_matched_through_boundary(
        sptr<colored_vertex_t> v1,
        sptr<colored_vertex_t> v2,
        std::string r,
        sptr<colored_vertex_t>* b1_p,
        sptr<colored_vertex_t>* b2_p,
        bool use_flagged_graph)
{
    typedef std::pair<sptr<colored_vertex_t>, sptr<colored_vertex_t>> cvp_t;
    static std::map<std::pair<cvp_t, std::string>, cvp_t> memo; // We'll memoize any data to keep this
                                        // fast for repeated requests.
    DecodingGraph::matrix_entry_t error_data;
    if (use_flagged_graph) {
        error_data = (*this)[r].get_error_chain_data_from_flagged_graph(v1, v2);
    } else {
        error_data = (*this)[r].get_error_chain_data(v1, v2);
    }
    if (!error_data.error_chain_runs_through_boundary)  return false;

    cvp_t v1_v2 = std::make_pair(v1, v2);
    cvp_t v2_v1 = std::make_pair(v2, v1);

    auto v1_v2_r = std::make_pair(v1_v2, r);
    if (memo.count(v1_v2_r)) {
        auto& boundaries = memo[v1_v2_r];
        *b1_p = boundaries.first;
        *b2_p = boundaries.second;
        return true;
    }
    // Otherwise, we have to compute the two boundaries.
    auto path = error_data.error_chain;
    // Ignore both endpoints.
    sptr<colored_vertex_t> b1 = nullptr;
    sptr<colored_vertex_t> b2 = nullptr;
    if (path[0] != v1) std::reverse(path.begin(), path.end());
    for (uint i = 1; i < path.size()-1; i++) {
        auto w = std::static_pointer_cast<colored_vertex_t>(path[i]);
        if (is_colored_boundary(w)) {
            if (b1 == nullptr) {
                b1 = w;
                b2 = w;
            } else {
                // Then, there are two different boundaries.
                b2 = w;
            }
        }
    }
    auto v2_v1_r = std::make_pair(v2_v1, r);
    memo[v1_v2_r] = std::make_pair(b1, b2);
    memo[v2_v1_r] = std::make_pair(b2, b1);
    
    *b1_p = b1;
    *b2_p = b2;
    return true;
}

std::set<sptr<colored_vertex_t>>
ColoredDecodingGraph::get_all_incident_vertices(const std::set<sptr<colored_edge_t>>& edge_list, std::string color) {
    std::set<sptr<colored_vertex_t>> incident;
    for (auto e : edge_list) {
        auto v1 = std::reinterpret_pointer_cast<colored_vertex_t>(e->src);
        auto v2 = std::reinterpret_pointer_cast<colored_vertex_t>(e->dst);
        if (v1->color == color) incident.insert(v1);
        if (v2->color == color) incident.insert(v2);
    }
    return incident;
}


ColoredDecodingGraph
to_colored_decoding_graph(const stim::Circuit& circuit, DecodingGraph::Mode mode) {
    ColoredDecodingGraph graph(mode);
    // We need to create three circuits, one for each color restriction.
    const std::string restrictions[] = {"rg", "rb", "gb"};

    std::array<stim::Circuit, 3> subcircuits;
    subcircuits.fill(stim::Circuit());

    // Furthermore, as we will then pass each subcircuit through the ErrorAnalyzer,
    // we must also track the detector mappings from the original circuit to the
    // subcircuits.
    std::array<uint64_t, 3> sub_detector_ctr;
    sub_detector_ctr.fill(0);
    uint64_t detector_ctr = 0;

    // labeled_det_t: sub detector id, restricted lattice (i.e. "rg")
    typedef std::pair<uint64_t, std::string> labeled_det_t;
    std::map<labeled_det_t, uint64_t> sub_detector_map;

    circuit.for_each_operation([&] (stim::Operation op) {
        std::string gate_name(op.gate->name);
        if (gate_name == "DETECTOR") {
            // Check the color of the detector.
            uint64_t detector = detector_ctr;
            int color_id = (int) op.target_data.args[0];

            std::string color;
            if (color_id >= 0) {
                color = int_to_color(color_id);
                for (uint i = 0; i < 3; i++) {
                    std::string r = restrictions[i];
                    if (r[0] == color[0] || r[1] == color[0]) {
                        // Then, this detector belongs to this restricted lattice.
                        labeled_det_t sub_detector = std::make_pair(sub_detector_ctr[i]++, r);
                        sub_detector_map[sub_detector] = detector;
                        // Add this operation to the corresponding subcircuit as well.
                        subcircuits[i].append_operation(op);
                    }
                }
            } else {
                color = "none";
            }
            // Finally, to avoid issues down the line, make the vertex for the
            // detector here.
            sptr<colored_vertex_t> v = std::make_shared<colored_vertex_t>();
            v->id = detector;
            v->color = color;
            graph.add_vertex(v);
            detector_ctr++;
        } else {
            for (auto& sc : subcircuits) sc.append_operation(op);
        }
    });

    // Now, update the DecodingGraph for each subcircuit.
    for (uint i = 0; i < 3; i++) {
        std::string r = restrictions[i];
        stim::Circuit& sc = subcircuits[i];

        stim::DetectorErrorModel dem = 
            stim::ErrorAnalyzer::circuit_to_detector_error_model(
                sc,
                true,   // decompose_errors
                true,   // fold loops
                false,  // allow gauge detectors
                1.0,    // approx disjoint errors threshold
                false,  // ignore decomposition failures
                false
            );
        // Nothing needs to be done for the detector callback as we have already
        // added the detector.
        detector_callback_t det_f = [] (uint64_t d, std::array<fp_t, N_COORD> coords) {};

        error_callback_t err_f =
            [&] (fp_t prob, std::vector<uint64_t> dets, std::set<uint> frames)
            {
                if (prob == 0 || dets.size() == 0 || dets.size() > 2) {
                    return;  // Zero error probability -- not an edge.
                }
                // First, convert the detectors (which are sub detectors right now)
                // to their true ids.
                for (auto& d : dets) {
                    labeled_det_t ld = std::make_pair(d, r);
                    d = sub_detector_map[ld];
                }

                if (dets.size() == 1) {
                    // We will be connecting to the boundary, so we must choose the correct one.
                    uint64_t d = dets[0];
                    auto v = graph.get_vertex(d);
                    std::string c = v->color;
                    if (r == "rg" || r == "gr") {
                        if (c == "r")   dets.push_back(GREEN_BOUNDARY_INDEX);
                        else            dets.push_back(RED_BOUNDARY_INDEX);
                    } else if (r == "rb" || r == "br") {
                        if (c == "r")   dets.push_back(BLUE_BOUNDARY_INDEX);
                        else            dets.push_back(RED_BOUNDARY_INDEX);
                    } else {
                        if (c == "g")   dets.push_back(BLUE_BOUNDARY_INDEX);
                        else            dets.push_back(GREEN_BOUNDARY_INDEX);
                    }
                }
                // Now, there are only two detectors.
                auto v1 = graph.get_vertex(dets[0]);
                auto v2 = graph.get_vertex(dets[1]);
                auto e = graph.get_edge(v1, v2);
                if (e != nullptr) {
                    fp_t old_prob = e->error_probability;
                    std::set<uint> old_frames = e->frames;
                    if (frames == old_frames) {
                        prob = prob * (1-old_prob) + old_prob * (1-prob);
                    }
                } else {
                    // Create new edge if it does not exist.
                    e = std::make_shared<colored_edge_t>();
                    e->src = std::static_pointer_cast<void>(v1);
                    e->dst = std::static_pointer_cast<void>(v2);
                    graph.add_edge(e);
                }
                fp_t edge_weight = (fp_t)log10((1-prob)/prob);
                e->edge_weight = edge_weight;
                e->error_probability = prob;
                e->frames = frames;
            };
        // Read the detector error model for the subcircuit.
        uint det_offset = 0;
        std::array<fp_t, N_COORD> coord_offset;
        coord_offset.fill(0);
        read_detector_error_model(dem, 1, det_offset, coord_offset, err_f, det_f);
    }
    // Analyze the original circuit to get observable flips.
    stim::DetectorErrorModel dem = 
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
            circuit,
            false,  // decompose_errors
            true,   // fold loops
            false,  // allow gauge detectors
            1.0,    // approx disjoint errors threshold
            true,   // ignore decomposition failures
            false
        );
    detector_callback_t det_f = [] (uint64_t d, std::array<fp_t, N_COORD> coords) {};

    error_callback_t err_f =
        [&] (fp_t prob, std::vector<uint64_t> dets, std::set<uint> frames)
        {
            if (dets.size() == 0) return;   // Nothing to do.
            if (dets.size() == 1) {
                // This vertex is adjacent to two boundaries of opposing color.
                auto v = graph.get_vertex(dets[0]);
                sptr<colored_vertex_t> b1;
                sptr<colored_vertex_t> b2;
                if (v->color == "r") {
                    b1 = graph.get_vertex(GREEN_BOUNDARY_INDEX);
                    b2 = graph.get_vertex(BLUE_BOUNDARY_INDEX);
                } else if (v->color == "g") {
                    b1 = graph.get_vertex(RED_BOUNDARY_INDEX);
                    b2 = graph.get_vertex(BLUE_BOUNDARY_INDEX);
                } else {
                    b1 = graph.get_vertex(RED_BOUNDARY_INDEX);
                    b2 = graph.get_vertex(GREEN_BOUNDARY_INDEX);
                }
                face_t fc = make_face(v, b1, b2);
                graph.face_frame_map[fc] = frames;
                graph.face_prob_map[fc] = prob;
            } else if (dets.size() == 2) {
                auto v1 = graph.get_vertex(dets[0]);
                auto v2 = graph.get_vertex(dets[1]);
                sptr<colored_vertex_t> b;
                if (v1->color == "r") {
                    if (v2->color == "g")   b = graph.get_vertex(BLUE_BOUNDARY_INDEX);
                    else                    b = graph.get_vertex(GREEN_BOUNDARY_INDEX);
                } else if (v1->color == "g") {
                    if (v2->color == "r")   b = graph.get_vertex(BLUE_BOUNDARY_INDEX);
                    else                    b = graph.get_vertex(RED_BOUNDARY_INDEX);
                } else {
                    if (v2->color == "r")   b = graph.get_vertex(GREEN_BOUNDARY_INDEX);
                    else                    b = graph.get_vertex(RED_BOUNDARY_INDEX);
                }
                face_t fc = make_face(v1, v2, b);
                graph.face_frame_map[fc] = frames;
                graph.face_prob_map[fc] = prob;
            } else {
                auto v1 = graph.get_vertex(dets[0]);
                auto v2 = graph.get_vertex(dets[1]);
                auto v3 = graph.get_vertex(dets[2]);
                face_t fc = make_face(v1, v2, v3);
                graph.face_frame_map[fc] = frames;
                graph.face_prob_map[fc] = prob;
            }
        };
    uint det_offset = 0;
    std::array<fp_t, N_COORD> coord_offset;
    coord_offset.fill(0);
    read_detector_error_model(dem, 1, det_offset, coord_offset, err_f, det_f);

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
                std::vector<uint64_t> detectors;
                std::set<uint> frames;
                 
                fp_t e_prob = (fp_t)inst.arg_data[0];
                for (stim::DemTarget target : inst.target_data) {
                    if (target.is_relative_detector_id()) {
                        // This is a detector, add it to the list.
                        detectors.push_back(
                                (uint64_t)target.data + det_offset);
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
