/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#include "qontra/graph/decoding_graph.h"

#include <vtils/utility.h>

#include <initializer_list>

namespace qontra {
namespace graph {

using namespace decoding;

DecodingGraph::DecodingGraph(const DetailedStimCircuit& circuit, size_t flips_per_error)
    :HyperGraph(),
    number_of_colors(circuit.number_of_colors_in_circuit),
    error_polynomial(),
    expected_errors(),
    dijkstra_graph_map(),
    flagged_dijkstra_graph_map(),
    distance_matrix_map(),
    flagged_distance_matrix_map(),
    active_flags(),
    flag_edge_map(),
    flag_detectors(circuit.flag_detectors),
    flags_are_active(false)
{
    stim::DetectorErrorModel dem =
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
                circuit,
                true,   // decompose errors
                true,   // fold loops
                false,  // allow gauge detectors (non-deterministic detectors)
                1.0,    // approx disjoint errors threshold
                true,   // ignore decomposition failures
                false
            );
    // Create boundaries first.
    if (number_of_colors == 0) {
        // Then we have a single boundary.
        uint64_t d = get_color_boundary_index(COLOR_ANY);
        sptr<vertex_t> vb = make_and_add_vertex(d);
        vb->is_boundary_vertex = true;
    } else {
        std::vector<sptr<vertex_t>> boundary_vertices;
        for (int c = 0; c < number_of_colors; c++) {
            uint64_t d = get_color_boundary_index(c);
            sptr<vertex_t> vb = make_and_add_vertex(d);
            vb->is_boundary_vertex = true;
            vb->color = c;
            boundary_vertices.push_back(vb);
        }
        // Make all of these boundaries connected to each other by a 1.0 probability hyperedge.
        sptr<hyperedge_t> e = make_and_add_edge(boundary_vertices);
        e->probability = 1.0;
    }
    auto df = 
        [&] (uint64_t d)
        {
            if (!this->contains(d)) this->make_and_add_vertex(d);
        };
    auto ef =
        [&] (fp_t p, std::vector<uint64_t> detectors, std::set<uint64_t> frames)
        {
            if (p == 0) return;
            // Remove all flags from the detectors.
            std::vector<uint64_t> flags;
            for (auto it = detectors.begin(); it != detectors.end(); ) {
                if (flag_detectors.count(*it)) {
                    flags.push_back(*it);
                    it = detectors.erase(it);
                } else {
                    it++;
                }
            }
            if (detectors.size() == 0) return;
            // If it is less than the expected number of flips, then add
            // boundaries.
            if (detectors.size() < flips_per_error) {
                if (number_of_colors == 0) {
                    uint64_t b = get_color_boundary_index(COLOR_ANY);
                    detectors.push_back(b);
                } else {
                    std::vector<int> colors_in_detectors;
                    for (uint64_t d : detectors) {
                        colors_in_detectors.push_back(circuit.detector_color_map.at(d));
                    }
                    std::vector<uint64_t> boundary_indices;
                    for (int c : get_complementary_colors_to(colors_in_detectors, number_of_colors)) {
                        boundary_indices.push_back(get_color_boundary_index(c));
                    }
                    // If boundary_indices.size() + dets.size() overflows, it's
                    // just best to leave dets as is (likely a measurement
                    // error).
                    if (boundary_indices.size() + detectors.size() == flips_per_error) {
                        vtils::push_back_range(detectors, boundary_indices);
                    }
                }
            }
            // If there is an overflow, print out debug info and exit.
            if (detectors.size() > flips_per_error) {
                std::cerr << "[ DecodingGraph ] found error flipping detectors [";
                for (uint64_t d : detectors) std::cerr << " " << d;
                std::cerr << " ] and flags [";
                for (uint64_t d : flags) std::cerr << " " << d;
                std::cerr << " ] which had more than " << flips_per_error << " flips." << std::endl;
                exit(1);
            }
            // Create hyperedge now.
            std::vector<sptr<vertex_t>> vlist;
            for (uint64_t d : detectors) {
                if (!this->contains(d)) {
                    vlist.push_back(make_and_add_vertex_(d, circuit));
                } else {
                    vlist.push_back(this->get_vertex(d));
                }
            }
            sptr<hyperedge_t> e = this->get_edge(vlist);
            if (e != nullptr) {
                fp_t r = e->probability;
                if (frames == e->frames) {
                    p = p*(1-r) + r*(1-p);
                }
            } else {
                e = this->make_edge(vlist);
            }
            e->probability = p;
            e->frames = frames;
            // If this is a flag edge, do not add it to the graph. Instead, add
            // it to the flag_edge_map.
            if (flags.empty()) {
                this->add_edge(e);
            } else {
                for (uint64_t f : flags) {
                    this->flag_edge_map[f].insert(e);
                }
            }
        };
    size_t detector_offset = 0;
    read_detector_error_model(dem, 1, detector_offset, ef, df);
}

void
DecodingGraph::dijkstra_(int c1, int c2, sptr<vertex_t> from) {
    if (c1 > c2) return dijkstra_(c2, c1, from);
    auto c1_c2 = std::make_pair(c1, c2);
    
    auto& dgr_map = flags_are_active ? flagged_dijkstra_graph_map : dijkstra_graph_map;
    auto& dm_map = flags_are_active ? flagged_distance_matrix_map : distance_matrix_map;

    if (!dgr_map.count(c1_c2)) {
        make_dijkstra_graph(c1, c2);
    }
    // Retrieve the dijkstra graph from the map. Place it back later.
    uptr<DijkstraGraph> dgr = std::move(dgr_map[c1_c2]);
    std::map<sptr<vertex_t>, fp_t> dist;
    std::map<sptr<vertex_t>, sptr<vertex_t>> pred;
    dijkstra(dgr.get(), from, dist, pred, 
        [] (sptr<edge_t> e)
        {
            return compute_weight(e->probability);
        });
    // Now, update the distance matrix. We'll do this manually to optimize for
    // speed.
    auto& dm = dm_map[c1_c2];
    for (sptr<vertex_t> w : dgr->get_vertices()) {
        error_chain_t ec;

        bool failed_to_converge = false;
        size_t n_iter = 0;

        sptr<vertex_t> curr = w;
        while (curr != from) {
            sptr<vertex_t> next = pred[curr];
            if (curr == next) {
                failed_to_converge = true;
                break;
            }
            sptr<edge_t> e = dgr->get_edge(curr, next);
            // Update error chain data.
            ec.length++;
            ec.probability *= e->probability;
            ec.path.push_back(curr);
            if (curr->is_boundary_vertex) {
                ec.runs_through_boundary = true;
                ec.boundary_vertices.push_back(curr);
            }
            curr = next;
            // Check if are taking too long.
            if ((++n_iter) >= 10'000) {
                failed_to_converge = true;
                break;
            }
        }
        // Some last updates (adding from).
        ec.weight = dist.at(w);
        ec.length++;
        ec.path.push_back(from);
        if (from->is_boundary_vertex) {
            ec.runs_through_boundary = true;
            ec.boundary_vertices.push_back(from);
        }
        // Handle the case where we fail to converge:
        if (failed_to_converge) {
            ec.length = std::numeric_limits<size_t>::max();
            ec.probability = 0.0;
            ec.weight = std::numeric_limits<fp_t>::max();
            ec.runs_through_boundary = false;
            ec.path.clear();
            ec.boundary_vertices.clear();
        }
        // Update dm.
        dm[from][w] = std::move(ec);
    }
    dgr_map[c1_c2] = std::move(dgr);
}

void
DecodingGraph::make_dijkstra_graph(int c1, int c2) {
    // We need to build a suitable graph for Dijkstra's:
    uptr<Graph<vertex_t, edge_t>> dgr = std::make_unique<Graph<vertex_t, edge_t>>();
    // Populate dgr.
    for (sptr<vertex_t> v : get_vertices()) {
        if (c1 == COLOR_ANY || v->color == c1 || v->color == c2) {
            dgr->add_vertex(v);
        }
    }
    // Now, add edges to the graph.
    // First handle flag edges.
    std::set<sptr<hyperedge_t>> visited_flag_edges;
    fp_t renorm_factor = 1.0;
    for (uint64_t fd : active_flags) {
        for (sptr<hyperedge_t> he : flag_edge_map[fd]) {
            if (visited_flag_edges.count(he)) continue;

            fp_t p = he->probability;
            for (size_t i = 0; i < he->get_order(); i++) {
                sptr<vertex_t> v = he->get<vertex_t>(i);
                if (!dgr->contains(v)) continue;
                for (size_t j = i+1; j < he->get_order(); j++) {
                    sptr<vertex_t> w = he->get<vertex_t>(j);
                    if (!dgr->contains(w)) continue;
                    // Make edge if it does not exist. Update probability.
                    sptr<edge_t> e = dgr->get_edge(v, w);
                    if (e == nullptr) {
                        e = dgr->make_and_add_edge(v, w);
                    }
                    fp_t& r = e->probability;
                    r = (1-r)*p + (1-p)*r;
                }
            }
            renorm_factor *= p;
            visited_flag_edges.insert(he);
        }
    }
    // Now handle other edges.
    for (sptr<vertex_t> v : dgr->get_vertices()) {
        for (sptr<vertex_t> w : dgr->get_vertices()) {
            if (v == w) continue;
            sptr<edge_t> e = dgr->get_edge(v, w);
            if (e == nullptr) {
                e = dgr->make_and_add_edge(v, w);
            }
            fp_t p = 0.0;
            for (sptr<hyperedge_t> he : get_common_hyperedges({v, w})) {
                fp_t r = he->probability;
                p = (1-p)*r + (1-r)*p;
            }
            e->probability = renorm_factor*p;
        }
    }
    // Now the graph is done.
    auto c1_c2 = std::make_pair(c1, c2);
    if (flags_are_active) {
        flagged_dijkstra_graph_map[c1_c2] = std::move(dgr);
    } else {
        dijkstra_graph_map[c1_c2] = std::move(dgr);
    }
}

void
DecodingGraph::build_error_polynomial() {
    sptr<hyperedge_t> e0 = edges[0];
    poly_t pX{1 - e0->probability, e0->probability};
    
    fp_t expectation = 0.0;
    for (size_t i = 1; i < edges.size(); i++) {
        sptr<hyperedge_t> e = edges[i];
        poly_t a(pX.size()+1, 0),
               b(pX.size()+1, 0);
        for (size_t j = 0; j < pX.size(); j++) {
            a[j] = pX[j] * (1 - e->probability);
            b[j] = pX[j] * e->probability;
        }
        pX = std::move(a);
        for (size_t j = 0; j < pX.size(); j++) {
            pX[j] += b[j];
            if (i == edges.size()-1) {
                expectation += j * pX[j];
            }
        }
    }
    error_polynomial = std::move(pX);
    expected_errors = expectation;
}

}   // graph
}   // qontra
