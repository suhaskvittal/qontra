/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#include "qontra/graph/decoding_graph.h"

#include <vtils/set_algebra.h>
#include <vtils/two_level_map.h>

#ifdef DECODER_PERF
#include <vtils/timer.h>
#endif

namespace qontra {
namespace graph {

using namespace decoding;

error_chain_t
DecodingGraph::get_error_chain(
        uint64_t d1,
        uint64_t d2,
        int c1,
        int c2,
        bool force_unflagged)
{
    return get_error_chain(get_vertex(d1), get_vertex(d2), c1, c2, force_unflagged);
}

error_chain_t
DecodingGraph::get_error_chain(
        sptr<decoding::vertex_t> v1,
        sptr<decoding::vertex_t> v2,
        int c1,
        int c2,
        bool force_unflagged) 
{
    update_state();
    if (c1 > c2) std::swap(c1, c2);
    auto c1_c2 = std::make_pair(c1, c2);
    auto& dm_map = flags_are_active && !force_unflagged 
                    ? flagged_distance_matrix_map[c1_c2] : distance_matrix_map[c1_c2];
    if (!dm_map.count(v1) || !dm_map.at(v1).count(v2)) {
        dijkstra_(c1, c2, v1);
    }
    return dm_map[v1][v2];
}

void
DecodingGraph::init_distances_for(int c1, int c2) {
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

fp_t
DecodingGraph::compute_renorm_factor(std::set<uint64_t> flags) {
    fp_t f = 1.0;

    std::set<uint64_t> remaining_flags = active_flags - flags;
    for (sptr<hyperedge_t> e : nod_edges) {
        if (remaining_flags.empty()) break;
        std::set<uint64_t> flag_intersect = e->flags * remaining_flags;
        const size_t nf = flag_intersect.size();
        const fp_t p = e->probability;
        if (nf == e->flags.size()) {
            remaining_flags -= flag_intersect;
            f *= p;
        }
    }
    return f;
}

fp_t
DecodingGraph::renormalized_edge_probability(sptr<hyperedge_t> e) {
    fp_t p = e->probability;
    // Three renormalizations.
    //
    // First: flag renormalization.
    if (e->flags.empty()) {
        p *= renorm_factor;
    } else {
        p *= compute_renorm_factor(e->flags);
        // Second (while we're still in this body), renormalize flags based
        // on order.
        p = pow(p, static_cast<fp_t>(e->get_order()-1));
    }
    if (!reweigh_for_detectors) return p;
    // Third: detector-based renormalization. Edges adjacent to an active
    // detector have higher probability.
    bool contains_boundary = false;
    size_t det_w = 0;
    for (sptr<vertex_t> v : e->get<vertex_t>()) {
        if (active_detectors.count(v->id)) {
            det_w++;
        }
        if (v->is_boundary_vertex && !contains_boundary) {
            // Check if a restricted lattice needs a boundary.
            for (int c1 = 0; c1 < number_of_colors; c1++) {
                for (int c2 = c1+1; c2 < number_of_colors; c2++) {
                    if (v->color != c1 && v->color != c2) continue;
                    size_t cnt = 0;
                    for (uint64_t d : active_detectors) {
                        if (d == get_color_boundary_index(COLOR_ANY)) continue;
                        sptr<vertex_t> w = get_vertex(d);
                        cnt += (w->color == c1 || w->color == c2);
                    }
                    if (cnt & 1) contains_boundary = true;
                }
            }
        }
    }
    if (det_w >= 2 || (det_w == 1 && contains_boundary)) {
        p = sqrt(p);
    } else if (det_w == 1) {
        p *= 10;
    }
    return p;
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
    std::cout << "[ DecodingGraph ] took " << t*1e-9 << "s to perform Dijkstra's for n = " 
        << dgr->n() << ", m = " << dgr->m() << std::endl;
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
    if (flags_are_active && !reweigh_for_detectors) {
        // Get the non-flagged map.
        auto& _dm = distance_matrix_map[c1_c2];
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
        // Now add the active detectors.
        for (uint64_t d : active_detectors) {
            if (d == get_color_boundary_index(COLOR_ANY)) {
                continue;
            } else {
                sptr<vertex_t> v = get_vertex(d);
                if (c1 == COLOR_ANY || v->color == c1 || v->color == c2) {
                    dgr->add_vertex(v);
                }
            }
        }
        // Also add all detectors in any flag edges.
        for (sptr<hyperedge_t> he : get_flag_edges()) {
            for (sptr<vertex_t> v : he->get<vertex_t>()) {
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
                    sptr<vertex_t> x = ec.path[i-1],
                                    y = ec.path[i];
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
    for (sptr<hyperedge_t> he : get_flag_edges()) {
        fp_t p = renormalized_edge_probability(he);
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
#ifdef MEMORY_DEBUG
                std::cout << "Added flag edge between " << print_v(v) << " and " << print_v(w) << ", P = " << r << std::endl;
#endif
            }
        }
    }
#ifdef MEMORY_DEBUG
    std::cout << "Renorm factor = " << renorm_factor << std::endl;
#endif
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
                e->probability = 0.5;
            } else {
                fp_t p = 0.0;
                for (sptr<hyperedge_t> he : get_common_hyperedges({v, w})) {
                    fp_t r = renormalized_edge_probability(he);
                    p = (1-p)*r + (1-r)*p;
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

}   // graph
}   // qontra
