/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#include "qontra/graph/decoding_graph.h"
#include "qontra/graph/decoding_graph/edge_class.h"

#include <vtils/set_algebra.h>
#include <vtils/utility.h>

#include <initializer_list>

#ifdef DECODER_PERF
#include <vtils/timer.h>
#endif

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
    base_probability_map(),
    active_flags(),
    flag_detectors(circuit.flag_detectors),
    edge_classes(),
    flag_class_map(),
    edge_class_map(),
    nod_edges(),
    all_edges(),
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
                /*
                std::cout << "Detectorless event: F[";
                for (uint64_t f : flags) std::cout << " " << f;
                std::cout << " ], frames =";
                for (uint64_t fr : frames) std::cout << " " << fr;
                std::cout << ", prob = " << p << std::endl;
                */

                sptr<hyperedge_t> e = this->make_edge({});
                e->flags = std::set<uint64_t>(flags.begin(), flags.end());
                e->frames = frames;
                e->probability = p;
                nod_edges.push_back(e);
                all_edges.push_back(e);
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

void
DecodingGraph::immediately_initialize_distances_for(int c1, int c2) {
#ifdef DECODER_PERF
    vtils::Timer timer;
    fp_t t;
#endif
    update_state();
    if (c1 > c2) std::swap(c1, c2);
    auto c1_c2 = std::make_pair(c1, c2);
    auto& dgr_map = flags_are_active ? flagged_dijkstra_graph_map : dijkstra_graph_map;
    auto& dm_map = flags_are_active ? flagged_distance_matrix_map : distance_matrix_map;

    if (!dgr_map.count(c1_c2)) {
        make_dijkstra_graph(c1, c2);
    }
    uptr<DijkstraGraph>& dgr = dgr_map[c1_c2];

    vtils::TwoLevelMap<sptr<vertex_t>, sptr<vertex_t>, fp_t> dist;
    vtils::TwoLevelMap<sptr<vertex_t>, sptr<vertex_t>, sptr<vertex_t>> pred;
#ifdef DECODER_PERF
    timer.clk_start();
#endif
    floyd_warshall(dgr.get(), dist, pred,
        [] (sptr<edge_t> e)
        {
            return compute_weight(e->probability);
        });
#ifdef DECODER_PERF
    t = timer.clk_end();
    std::cout << "[ DecodingGraph ] took " << t*1e-9 << "s to run Floyd-Warshall for n = "
        << dgr->n() << " and m = " << dgr->m() << std::endl;
#endif
    auto& dm = dm_map[c1_c2];
    for (sptr<vertex_t> v : dgr->get_vertices()) {
        update_paths(dgr, dm, v, dgr->get_vertices(), dist[v], pred[v]);
    }
}

sptr<hyperedge_t>
DecodingGraph::get_best_shared_edge(std::vector<sptr<vertex_t>> vlist) {
    auto common = get_common_hyperedges(vlist);
    if (common.empty()) return nullptr;
    if (common.size() == 1) return common[0];
    // Choose hyperedge with lowest order. Break ties with probability.
    auto it = std::min_element(common.begin(), common.end(),
                [] (sptr<hyperedge_t> x, sptr<hyperedge_t> y)
                {
                    return (x->get_order() < y->get_order())
                        || (x->get_order() == y->get_order() && x->probability > y->probability);
                });
    return *it;
}

sptr<hyperedge_t>
DecodingGraph::get_base_edge(sptr<hyperedge_t> e) {
    if (!edge_class_map.count(e)) {
        return e;
    }
    // We want to (1) get the representative of e's class, (2) flatten the rep,
    // and (3) get the class of the flattened rep. Then, we will search the class
    // for an edge with identical frames to e and an identical number of flags to e.
    const EdgeClass& c = edge_class_map.at(e);
    sptr<hyperedge_t> rep = c.get_representative();
    // Cannot flatten an edge represented by a flag edge.
    if (rep->flags.size()) {
        return e;
    }

    std::vector<sptr<vertex_t>> vlist;
    for (size_t i = 0; i < rep->get_order(); i++) {
        vlist.push_back(rep->get<vertex_t>(i)->get_base());
    }
    sptr<hyperedge_t> frep = get_edge(vlist);
    const EdgeClass& fc = edge_class_map.at(frep);
    for (sptr<hyperedge_t> fe : fc.get_edges()) {
        if (fe->flags.size() == e->flags.size() && fe->frames == e->frames) {
            return fe;
        }
    }
    // Could not flatten this.
    return e;
}

std::vector<sptr<hyperedge_t>>
DecodingGraph::get_flag_edges() {
    // Select one flag edge from each equivalence class that has the most in common
    // with the active flags.
    std::vector<sptr<hyperedge_t>> edge_list;
    for (EdgeClass& c : edge_classes) {
        sptr<hyperedge_t> best_edge = get_best_flag_edge(c.get_edges());
        if (best_edge != nullptr) {
            edge_list.push_back(best_edge);
        }
    }
    return edge_list;
}

std::map<sptr<hyperedge_t>, sptr<hyperedge_t>>
DecodingGraph::get_best_rep_map() {
    std::map<sptr<hyperedge_t>, sptr<hyperedge_t>> rep_map;
    for (uint64_t f : active_flags) {
        for (const EdgeClass& c : flag_class_map[f]) {
            sptr<hyperedge_t> best_e = get_best_flag_edge(c.get_edges());
            if (best_e == nullptr) {
                best_e = c.get_representative();
            }
            rep_map[c.get_representative()] = best_e;
        }
    }
    return rep_map;
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
                if (number_of_colors == 0) {
                    if (rep->get_order() + 1 == flips_per_error) {
                        sptr<vertex_t> vb = get_boundary_vertex(COLOR_ANY);
                        c.add_vertex(vb);
                    }
                } else {
                    std::vector<sptr<vertex_t>> boundary_list = 
                        get_complementary_boundaries_to(rep->get<vertex_t>());
                    if (boundary_list.size() + rep->get_order() == flips_per_error) {
                        for (sptr<vertex_t> vb : boundary_list) c.add_vertex(vb);
                    }
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

        /*
        std::cout << "\nEquivalence class of D[";
        for (size_t i = 0; i < rep->get_order(); i++) {
            std::cout << " " << print_v(rep->get<vertex_t>(i));
        }
        std::cout << " ], F[";
        for (uint64_t f : rep->flags) {
            std::cout << " " << f;
        }
        std::cout << " ]" << std::endl;
        */

        // Update any edge probabilities.
        std::set<uint64_t> all_flags;
        for (sptr<hyperedge_t> e : c.get_edges()) {
            if (e->flags.size() && rep->flags.empty()) {
                // Increase the probability of this edge, as it is also a non flag edge.
                e->probability += rep->probability;
            } else if (e->flags.empty()) {
                // Check if a similar edge already exists. If so, update that.
                sptr<hyperedge_t> _e = get_edge(e->get<vertex_t>());
                if (_e != nullptr) {
                    fp_t& p1 = e->probability,
                        & p2 = _e->probability;
                    if (e->frames == _e->frames) {
                        p2 = (1-p1)*p2 + (1-p2)*p1;
                    }
                } else {
                    add_edge(e);
                }
            }
            vtils::insert_range(all_flags, e->flags);
            edge_class_map[e] = c;

            /*
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
            */
            all_edges.push_back(e);
        }
        edge_classes.push_back(c);
        for (uint64_t f : all_flags) {
            flag_class_map[f].push_back(c);
        }
    }
}

sptr<hyperedge_t>
DecodingGraph::get_best_flag_edge(std::vector<sptr<hyperedge_t>> edge_list) {
    sptr<hyperedge_t> best_edge = nullptr;
    size_t best_nf = 0;
    fp_t best_p = 0.0;

    for (sptr<hyperedge_t> e : edge_list) {
        std::set<uint64_t> flag_intersect = e->flags * active_flags;
        const size_t nf = flag_intersect.size();
        const fp_t p = e->probability;
        if (nf == e->flags.size() && 
            (nf > best_nf || (nf > 0 && nf == best_nf && p > best_p)))
        {
            best_edge = e;
            best_nf = nf;
            best_p = p;
        }
    }
    return best_edge;
}

void
DecodingGraph::dijkstra_(int c1, int c2, sptr<vertex_t> from, sptr<vertex_t> to) {
#ifdef DECODER_PERF
    vtils::Timer timer;
    fp_t t;
#endif

    if (c1 > c2) return dijkstra_(c2, c1, from, to);
    auto c1_c2 = std::make_pair(c1, c2);
    
    auto& dgr_map = flags_are_active ? flagged_dijkstra_graph_map : dijkstra_graph_map;
    auto& dm_map = flags_are_active ? flagged_distance_matrix_map : distance_matrix_map;

    if (!dgr_map.count(c1_c2)) {
        make_dijkstra_graph(c1, c2);
    }
    // Retrieve the dijkstra graph from the map. Place it back later.
    uptr<DijkstraGraph>& dgr = dgr_map[c1_c2];
    std::map<sptr<vertex_t>, fp_t> dist;
    std::map<sptr<vertex_t>, sptr<vertex_t>> pred;
#ifdef DECODER_PERF
    timer.clk_start();
#endif
    dijkstra(dgr.get(), from, dist, pred,
        [] (sptr<edge_t> e)
        {
            return compute_weight(e->probability);
        }, to);
#ifdef DECODER_PERF
    t = timer.clk_end();
    std::cout << "[ DecodingGraph ] took " << t*1e-9 << "s to perform Dijkstra's" << std::endl;
#endif
    // Now, update the distance matrix. We'll do this manually to optimize for
    // speed.
    auto& dm = dm_map[c1_c2];
    std::vector<sptr<vertex_t>> to_list = //dgr->get_vertices();
        to == nullptr ? dgr->get_vertices() : std::vector<sptr<vertex_t>>{to};
    update_paths(dgr, dm, from, to_list, dist, pred);
}

void
DecodingGraph::make_dijkstra_graph(int c1, int c2) {
    auto c1_c2 = std::make_pair(c1, c2);
    auto& dgr_map = flags_are_active ? flagged_dijkstra_graph_map : dijkstra_graph_map;

    // We need to build a suitable graph for Dijkstra's:
    uptr<Graph<vertex_t, edge_t>> dgr = std::make_unique<Graph<vertex_t, edge_t>>();
    // Populate dgr. If flags_are_active, then assume dgr is transient only add
    // vertices relevant to the active detectors.
    if (flags_are_active) {
        // Get the non-flagged map.
        auto& _dm = distance_matrix_map[c1_c2];
        // Now add the active detectors.
        for (uint64_t d : active_detectors) {
            if (d == get_color_boundary_index(COLOR_ANY)) {
                // Add the boundaries of interest.
                if (c1 == COLOR_ANY) {
                    sptr<vertex_t> vb = get_boundary_vertex(COLOR_ANY);
                    dgr->add_vertex(vb);
                } else {
                    sptr<vertex_t> vb1 = get_boundary_vertex(c1),
                                    vb2 = get_boundary_vertex(c2);
                    dgr->add_vertex(vb1);
                    dgr->add_vertex(vb2);
                }
            } else {
                sptr<vertex_t> v = get_vertex(d);
                if (c1 == COLOR_ANY || v->color == c1 || v->color == c2) {
                    dgr->add_vertex(v);
                }
            }
        }
        // Finally add all detectors in a path between any two vertices. Also
        // add edges preemptively.
        for (sptr<vertex_t> v : dgr->get_vertices()) {
            for (sptr<vertex_t> w : dgr->get_vertices()) {
                // Get the path in the non-flagged map. Add all detectors in the
                // path between v and w.
                error_chain_t ec = _dm[v][w];
                if (ec.path.empty()) continue;
                for (size_t i = 1; i < ec.path.size(); i++) {
                    sptr<vertex_t> x = ec.path[i-1];
                    sptr<vertex_t> y = ec.path[i-1];
                    if (!dgr->contains(y)) dgr->add_vertex(y);
                    dgr->make_and_add_edge(x, y);
                }
            }
        }
    } else {
        for (sptr<vertex_t> v : get_vertices()) {
            if (c1 == COLOR_ANY || v->color == c1 || v->color == c2) {
                dgr->add_vertex(v);
            }
        }
    }
    // Now, add edges to the graph.
    // First handle flag edges.
#ifdef DECODER_PERF
    vtils::Timer timer;
    fp_t t;

    timer.clk_start();
#endif
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
#ifdef DECODER_PERF
    t = timer.clk_end();
    std::cout << "[ DecodingGraph ] took " << t*1e-9 << "s to add flag edges to Dijkstra graph" << std::endl;
#endif

#ifdef DECODER_PERF
    timer.clk_start();
#endif
    // Now handle other edges.
    auto vertices = dgr->get_vertices();
    for (size_t i = 0; i < vertices.size(); i++) {
        sptr<vertex_t> v = vertices.at(i);
        for (size_t j = i+1; j < vertices.size(); j++) {
            sptr<vertex_t> w = vertices.at(j);
            sptr<edge_t> e = dgr->get_edge(v, w);
            if (e == nullptr) {
                e = dgr->make_and_add_edge(v, w);
            }
            if (v->is_boundary_vertex && w->is_boundary_vertex) {
                e->probability = 1.0;
            } else {
                fp_t p = 0.0;
                if (base_probability_map.count(v) && base_probability_map.at(v).count(w)) {
                    p = base_probability_map.at(v).at(w);
                } else {
                    for (sptr<hyperedge_t> he : get_common_hyperedges({v, w})) {
                        fp_t r = he->probability;
                        p = (1-p)*r + (1-r)*p;
                    }
                    base_probability_map[v][w] = p;
                    base_probability_map[w][v] = p;
                }
                p *= renorm_factor;
                e->probability = (1-e->probability)*p + (1-p)*e->probability;
            }
        }
    }
    // Finally, delete all vertices with degree 0.
    for (sptr<vertex_t> v : dgr->get_vertices()) {
        if (dgr->get_degree(v) == 0) dgr->delete_vertex(v);
    }
#ifdef DECODER_PERF
    t = timer.clk_end();
    std::cout << "[ DecodingGraph ] took " << t*1e-9 << "s to add normal edges to Dijkstra graph" << std::endl;
#endif
    // Now the graph is done.
    dgr_map[c1_c2] = std::move(dgr);
}

void
DecodingGraph::update_paths(
        uptr<DijkstraGraph>& dgr,
        DistanceMatrix<vertex_t, error_chain_t>& dm,
        sptr<vertex_t> from,
        std::vector<sptr<vertex_t>> to_list,
        const std::map<sptr<vertex_t>, fp_t>& dist,
        const std::map<sptr<vertex_t>, sptr<vertex_t>>& pred)
{
    for (sptr<vertex_t> w : to_list) {
        error_chain_t ec;

        bool failed_to_converge = false;
        size_t n_iter = 0;

        sptr<vertex_t> curr = w;
        while (curr != from) {
            sptr<vertex_t> next = pred.at(curr);
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
