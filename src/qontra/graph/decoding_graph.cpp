/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#include "qontra/graph/decoding_graph.h"
#include "qontra/graph/decoding_graph/edge_class.h"

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
                false,  // decompose errors
                true,   // fold loops
                false,  // allow gauge detectors (non-deterministic detectors)
                0.0,    // approx disjoint errors threshold
                false,  // ignore decomposition failures
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
    // Do not add edges to the graph immediately. First, we will collect all of them, and then
    // analyze the edges and add them accordingly (see resolve_edges).
    std::vector<sptr<hyperedge_t>> tentative_edges;
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
            if (detectors.size() == 0) {
                std::cout << "Detectorless event: F[";
                for (uint64_t f : flags) std::cout << " " << f;
                std::cout << " ], frames =";
                for (uint64_t fr : frames) std::cout << " " << fr;
                std::cout << ", prob = " << p << std::endl;
                return;
            }

            std::vector<sptr<vertex_t>> vlist;
            for (uint64_t d : detectors) {
                if (!this->contains(d)) {
                    vlist.push_back(make_and_add_vertex_(d, circuit));
                } else {
                    vlist.push_back(this->get_vertex(d));
                }
            }
            sptr<hyperedge_t> e = make_edge(vlist);
            e->probability = p;
            e->frames = frames;
            e->flags = std::set<uint64_t>(flags.begin(), flags.end());
            tentative_edges.push_back(e);
        };
    size_t detector_offset = 0;
    read_detector_error_model(dem, 1, detector_offset, ef, df);

    resolve_edges(tentative_edges, flips_per_error);
}

sptr<vertex_t>
DecodingGraph::make_and_add_vertex_(uint64_t d, const DetailedStimCircuit& circuit) {
    sptr<vertex_t> x = make_and_add_vertex(d);
    if (circuit.detector_color_map.count(d)) {
        x->color = circuit.detector_color_map.at(d);
    }
    if (circuit.detector_base_map.count(d)) {
        uint64_t bd = circuit.detector_base_map.at(d);
        if (d != bd && !this->contains(bd)) {
            sptr<vertex_t> _x = make_and_add_vertex_(bd, circuit);
            x->base = _x;
        } else {
            x->base = this->get_vertex(bd);
        }
    }
    return x;
}

void
DecodingGraph::resolve_edges(const std::vector<sptr<hyperedge_t>>& edge_list, size_t flips_per_error) {
    std::vector<EdgeClass> classes = EdgeClass::from_edges(edge_list);
    // Go through and add boundaries to any edges if necessary.
    for (EdgeClass& c : classes) {
        // Add boundaries based on whether or not the representative is a flag edge.
        sptr<hyperedge_t> rep = c.get_representative();
        if (rep->flags.empty()) {
            if (rep->get_order() < flips_per_error) { 
                std::vector<sptr<vertex_t>> boundary_list = 
                    get_complementary_boundaries_to(rep->get<vertex_t>());
                if (boundary_list.size() + rep->get_order() == flips_per_error) {
                    for (sptr<vertex_t> vb : boundary_list) c.add_vertex(vb);
                }
            }
        } else {
            // Here, we operate under the assumption that the representative edge is a two qubit error
            // that will flip two detectors.
            if (rep->get_order() == 1) {
                // We have a single detector that must be linked to a boundary. We can just link to the
                // boundary of the same color.
                int color = rep->get<vertex_t>(0)->color;
                c.add_vertex(get_boundary_vertex(color));
            }
        }

        std::cout << "\nEquivalence class of D[";
        for (size_t i = 0; i < rep->get_order(); i++) {
            std::cout << " " << print_v(rep->get<vertex_t>(i));
        }
        std::cout << " ], F[";
        for (uint64_t f : rep->flags) {
            std::cout << " " << f;
        }
        std::cout << " ]" << std::endl;

        for (sptr<hyperedge_t> e : c.get_edges()) {
            if (e->flags.size() && rep->flags.empty()) {
                // Increase the probability of this edge, as it is also a non flag edge.
                e->probability += rep->probability;
            }

            std::cout << "\tD[";
            for (size_t i = 0; i < e->get_order(); i++) {
                std::cout << " " << print_v(e->get<vertex_t>(i));
            }
            std::cout << " ], F[";
            for (uint64_t f : e->flags) {
                std::cout << " " << f;
            }
            std::cout << " ], frames:";
            for (uint64_t fr : e->frames) {
                std::cout << " " << fr;
            }
            std::cout << ", prob = " << e->probability <<  std::endl;

            safe_add_edge(e);
        }
    }
}

void
DecodingGraph::safe_add_edge(sptr<hyperedge_t> e) {
    if (e->flags.empty()) {
        // Then check if the edge already exists in the graph.
        sptr<hyperedge_t> _e = get_edge(e->get<vertex_t>());
        if (_e == nullptr) {
            add_edge(e);
        } else {
            if (e->frames == _e->frames) {
                fp_t& p1 = e->probability,
                    & p2 = _e->probability;
                p2 = (1-p1)*p2 + (1-p2)*p1;
            }
        }
    } else {
        for (uint64_t f : e->flags) {
            flag_edge_map[f].insert(e);
        }
    }
}

std::vector<sptr<hyperedge_t>>
DecodingGraph::get_flag_edges() {
    std::map<sptr<hyperedge_t>, std::set<uint64_t>> remaining_flags_map;
    // edge_list should only contain edges such that all flags associated with
    // the flag edge are present in the flag detectors.
    std::vector<sptr<hyperedge_t>> edge_list;
    for (uint64_t fd : active_flags) {
        for (sptr<hyperedge_t> he : flag_edge_map.at(fd)) {
            // We will only consider he once all of its flags are found.
            if (!remaining_flags_map.count(he)) {
                remaining_flags_map[he] = he->flags;
            }
            remaining_flags_map[he].erase(fd);
            if (remaining_flags_map[he].empty()) {
                edge_list.push_back(he);
            }
        }
    }
    return edge_list;
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
    fp_t renorm_factor = 1.0;
    for (sptr<hyperedge_t> he : get_flag_edges()) {
        fp_t p = he->probability;
        bool any_flag_edge = false;
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
                any_flag_edge = true;
            }
        }
        if (any_flag_edge) {
            renorm_factor *= p;
        }
    }

    std::cout << "renorm factor: " << renorm_factor << std::endl;

    // Now handle other edges.
    for (sptr<vertex_t> v : dgr->get_vertices()) {
        for (sptr<vertex_t> w : dgr->get_vertices()) {
            if (v <= w) continue;
            sptr<edge_t> e = dgr->get_edge(v, w);
            if (e == nullptr) {
                e = dgr->make_and_add_edge(v, w);
            }
            if (v->is_boundary_vertex && w->is_boundary_vertex) {
                e->probability = 1.0;
            } else {
                fp_t p = 0.0;
                for (sptr<hyperedge_t> he : get_common_hyperedges({v, w})) {
                    fp_t r = he->probability;
                    p = (1-p)*r + (1-r)*p;
                }
                p *= renorm_factor;
                e->probability = (1-e->probability)*p + (1-p)*e->probability;
            }
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
