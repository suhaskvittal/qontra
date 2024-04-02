/* 
 *  author: Suhas Vittal
 *  date:   17 February 2024
 * */

#include "qontra/decoder/restriction.h"

#include <vtils/set_algebra.h>
#include <vtils/utility.h>

#include <initializer_list>

const int COLOR_RED = 0;
const int COLOR_GREEN = 1;

namespace qontra {

using namespace graph;
using namespace decoding;

inline void
push_back_assignment(std::vector<assign_t>& arr, const assign_t& m) {
    if (m.w == nullptr || m.v == m.w) {
        return;
    }
    arr.push_back(m);
}

inline void
erase_from_incidence_map(const vpair_t& e, std::map<vpair_t, size_t>& incidence_map) {
    if (incidence_map.count(e) && (--incidence_map[e]) == 0) {
        incidence_map.erase(e);
    }
}

template <class SETLIKE> inline size_t
locally_matches(SETLIKE s1, const std::set<vpair_t>& s2, sptr<vertex_t> v) {
    size_t m = 0;
    for (const vpair_t& e : s2) {
        if (!s1.count(e)) return 0;
        else              m++;
    }
    return m;
}

inline bool
update_best_boundary(
        int intersection,
        fp_t& best_log_p,
        std::set<vpair_t>& best_boundary,
        stim::simd_bits_range_ref<SIMD_WIDTH> best_corr,
        fp_t log_p,
        const std::set<vpair_t>& boundary,
        stim::simd_bits_range_ref<SIMD_WIDTH> corr)
{
    if (intersection > 0 && log_p > best_log_p) {
        best_log_p = log_p;
        best_boundary = boundary;
        best_corr.clear();
        best_corr |= corr;
        return true;
    }
    return false;
}

inline void
update_correction(
        std::map<vpair_t, size_t>& incidence_map,
        stim::simd_bits_range_ref<SIMD_WIDTH> corr,
        fp_t& corr_log_pr,
        stim::simd_bits_range_ref<SIMD_WIDTH> local_corr,
        const std::set<vpair_t>& local_boundary,
        fp_t local_log_pr)
{
    corr ^= local_corr;
    for (const vpair_t& e : local_boundary) {
        erase_from_incidence_map(e, incidence_map);
    }
    corr_log_pr += local_log_pr;
}

bool
remove_widowed_edges(std::map<vpair_t, size_t>& incidence_map) {
    std::map<sptr<vertex_t>, size_t> vertex_inc_map;
    for (const auto& [e, cnt] : incidence_map) {
        const auto& [v1, v2] = e;
        vertex_inc_map[v1]++;
        vertex_inc_map[v2]++;
    }
    // Remove any pairs of vertices where both endpoints only have a single
    // incidence.
    bool any_removed = false;
    for (auto it = incidence_map.begin(); it != incidence_map.end(); ) {
        const auto& [v1, v2] = it->first;
        if (vertex_inc_map.at(v1) == 1 && vertex_inc_map.at(v2) == 1) {
#ifdef MEMORY_DEBUG
            std::cout << "removed widowed edge [ " << print_v(v1) << " " << print_v(v2) << " ], cnt = "
                << it->second << std::endl;
#endif
            any_removed |= (it->second % 2 == 1);
            it = incidence_map.erase(it);
        } else {
            it++;
        }
    }
    return any_removed;
}

Decoder::result_t
RestrictionDecoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    // Reset data structures:
    in_cc_map.clear();
    not_cc_map.clear();
    component_edge_sets.clear();
    triggered_flag_edges.clear();

    const size_t n_obs = circuit.count_observables();
    stim::simd_bits<SIMD_WIDTH> corr(n_obs);
    load_syndrome(syndrome);
#ifdef MEMORY_DEBUG
    std::cout << "syndrome: D[";
    for (uint64_t d : detectors) std::cout << " " << d;
    std::cout << " ], F[";
    for (uint64_t f : flags) std::cout << " " << f;
    std::cout << " ]" << std::endl;
#endif

    if (detectors.empty()) return ret_no_detectors();

    corr ^= get_base_corr();

    // Compute the MWPM for each restricted lattice.
#ifdef DECODER_PERF
    timer.clk_start();
    fp_t t;
#endif
    std::map<sptr<hyperedge_t>, size_t> flag_edge_ctr_map;
    auto matchings = compute_matchings(syndrome);
    // Update incidence map.
    for (const assign_t& m : matchings) {
        for (sptr<hyperedge_t> e : m.flag_edges) {
            if ((++flag_edge_ctr_map[e]) == 2) {
#ifdef MEMORY_DEBUG
                std::cout << "Applying flag edge [";
                for (sptr<vertex_t> v : e->get<vertex_t>()) {
                    std::cout << " " << print_v(v);
                }
                std::cout << " ]" << std::endl;
#endif
                for (uint64_t fr : e->frames) corr[fr] ^= 1;
            }
        }
    }
    // After retrieving the matchings, we need to resolve any flag edges. If a flag edge appears in
    // two restricted lattices, immediately apply its corrections and remove it from any paths.
    //
    // Split any boundary matchings as well.
    std::vector<assign_t> new_matchings;
    for (const assign_t& m : matchings) {
        split_assignment(new_matchings, m, flag_edge_ctr_map);
    }
    matchings = std::move(new_matchings);
#ifdef MEMORY_DEBUG
    for (assign_t x : matchings) {
        std::cout << "L(" << x.c1 << "," << x.c2 << "):" << print_v(x.v) << 
            " <---> " << print_v(x.w) << ", path:";
        for (sptr<vertex_t> v : x.path) std::cout << " " << print_v(v);
        std::cout << std::endl;
    }
#endif
#ifdef DECODER_PERF
    t = timer.clk_end();
    std::cout << "[ RestrictionDecoder ] took " << t*1e-9 << "s to match restricted lattices" << std::endl;
#endif
    // Compute connected components.
    std::vector<component_t> components = compute_connected_components(matchings);
    // Identify all edges inside a connected component (and which component it is), and
    // do the same for all edges outside a connected component. Furthermore, maintain a
    // counter for the number of times an edge appears in a connected component.
    // 
    // We need to track all matches that are in a component to initialize not_cc_map.
    std::set<sptr<hyperedge_t>> applied_flag_edges;

    std::set<assign_t> in_cc_assignments;
    for (const component_t& cc : components) {
        int color = cc.color;
        std::set<vpair_t> edge_set;
        for (const assign_t& m : cc.assignments) {
            for (sptr<hyperedge_t> e : m.flag_edges) {
                if (applied_flag_edges.count(e)) continue;
                for (uint64_t fr : e->frames) {
                    corr[fr] ^= 1;
                }
                applied_flag_edges.insert(e);
#ifdef MEMORY_DEBUG
                std::cout << "Applied flag edge [";
                for (sptr<vertex_t> v : e->get<vertex_t>()) {
                    std::cout << " " << print_v(v);
                }
                std::cout << " ], frames =";
                for (uint64_t fr : e->frames) std::cout << " " << fr;
                std::cout << std::endl;
#endif
            }
            in_cc_assignments.insert(m);
            std::set<vpair_t> tmp = insert_error_chain_into(in_cc_map, m.path, color, m.c1, m.c2);
            vtils::insert_range(edge_set, tmp);
        }
        component_edge_sets.emplace_back(edge_set, color);
    }
    // Do not_cc now.
    for (const assign_t& m : matchings) {
        if (in_cc_assignments.count(m)) continue;
        if (m.c1 != COLOR_RED && m.c2 != COLOR_RED) continue;
        for (sptr<hyperedge_t> e : m.flag_edges) {
            if (applied_flag_edges.count(e)) continue;
            for (uint64_t fr : e->frames) {
                corr[fr] ^= 1;
            }
            applied_flag_edges.insert(e);
#ifdef MEMORY_DEBUG
            std::cout << "Applied flag edge [";
            for (sptr<vertex_t> v : e->get<vertex_t>()) {
                std::cout << " " << print_v(v);
            }
            std::cout << " ], frames =";
            for (uint64_t fr : e->frames) std::cout << " " << fr;
            std::cout << std::endl;
#endif
        }
        if (m.v == nullptr) {
            continue;
        }
        insert_error_chain_into(not_cc_map, m.path, COLOR_RED, m.c1, m.c2);
    }

    if (in_cc_map.empty() && not_cc_map.empty()) {
        return { 0, corr };
    }

    // Compute the best_rep_map for use in lifting (pass by ref).
    auto best_rep_map = decoding_graph->get_best_rep_map();
    best_rep_map = flatten_edge_map(best_rep_map);
    // Perform the lifting procedure. Twice if we have any triggered flag edges.
    // Copy incidence data structures as these will be modified by lifting.
    std::map<vpair_t, size_t> _in_cc_map(in_cc_map),
                              _not_cc_map(not_cc_map);
    std::vector<std::pair<std::set<vpair_t>, int>> _c_edge_sets(component_edge_sets);

    stim::simd_bits<SIMD_WIDTH> corr1(corr),
                                corr2(corr);
    fp_t log_p1 = lifting(corr1, best_rep_map);
#ifdef MEMORY_DEBUG
    std::cout << "corr1 = ";
    for (size_t i = 0; i < n_obs; i++) std::cout << corr1[i]+0;
    std::cout << std::endl;
#endif
    // If there are no triggered flag edges. Return here.
    if (triggered_flag_edges.empty()) {
        return { 0.0, corr1 };
    }
#ifdef MEMORY_DEBUG
    std::cout << "Triggered flag edges:" << std::endl;
    for (auto& [e, path, map_ref] : triggered_flag_edges) {
        std::cout << "\t[";
        for (sptr<vertex_t> v : e->get<vertex_t>()) std::cout << " " << print_v(v->get_base());
        std::cout << " ], frames =";
        for (uint64_t fr : e->frames) std::cout << " " << fr;
        std::cout << ", path = [";
        for (sptr<vertex_t> v : path) std::cout << " " << print_v(v->get_base());
        std::cout << " ]" << std::endl;
    }
#endif
    // Otherwise, perform the lifting procedure again, this time removing all edges corresponding
    // to triggered flag edges. Also update corr for each removed flag edge.
    in_cc_map = std::move(_in_cc_map);
    not_cc_map = std::move(_not_cc_map);
    component_edge_sets = std::move(_c_edge_sets);
    fp_t log_p2 = 0.0;

    std::set<sptr<hyperedge_t>> visited_flag_edges;
    for (const auto& [he, vlist, inc_map_ref] : triggered_flag_edges) {
        for (size_t i = 1; i < vlist.size(); i++) {
            sptr<vertex_t> v = vlist.at(i-1)->get_base(),
                           w = vlist.at(i)->get_base();
            vpair_t e = make_vpair(v, w);
            erase_from_incidence_map(e, inc_map_ref);
        }
        // Check if there is a similar hyperedge.
        bool found_similar = false;
        // Apply the edge's frame changes.
        if (!found_similar) {
            log_p2 += log(he->probability);
            for (uint64_t fr : he->frames) corr2[fr] ^= 1;
            visited_flag_edges.insert(he);
        }
    }
    log_p2 += lifting(corr2, best_rep_map);
#ifdef MEMORY_DEBUG
    std::cout << "corr1 = ";
    for (size_t i = 0; i < n_obs; i++) std::cout << corr1[i]+0;
    std::cout << std::endl;

    std::cout << "corr2 = ";
    for (size_t i = 0; i < n_obs; i++) std::cout << corr2[i]+0;
    std::cout << std::endl;

    std::cout << "log probs: " << log_p1 << " , " << log_p2 << std::endl;
    if (corr1 != corr2) std::cout << "correction mismatch detected.\n";
#endif
    corr = std::move(log_p1 > log_p2 ? corr1 : corr2);
    return { 0.0, corr };
}

std::vector<assign_t>
RestrictionDecoder::compute_matchings(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    std::vector<assign_t> matchings;
    for (int c1 = 0; c1 < decoding_graph->number_of_colors; c1++) {
        for (int c2 = c1+1; c2 < decoding_graph->number_of_colors; c2++) {
#ifdef MEMORY_DEBUG
            std::cout << "Matchings on L(" << c1 << ", " << c2 << ")-----------" << std::endl;
#endif
            load_syndrome(syndrome, c1, c2, false);
            std::vector<assign_t> _matchings = compute_matching(c1, c2);
            vtils::push_back_range(matchings, _matchings);
#ifdef MEMORY_DEBUG
            for (assign_t x : _matchings) {
                std::cout << "\t" << print_v(x.v) << " <---> " << print_v(x.w) << ", path:";
                for (sptr<vertex_t> v : x.path) std::cout << " " << print_v(v);
                std::cout << std::endl;
            }
#endif
        }
    }
    return matchings;
}

void
RestrictionDecoder::split_assignment(
        std::vector<assign_t>& assign_arr,
        const assign_t& m,
        const std::map<sptr<hyperedge_t>, size_t>& flag_ctr_map)
{
    assign_t curr;
    curr.c1 = m.c1;
    curr.c2 = m.c2;
    curr.v = m.v;
    curr.w = nullptr;
    curr.path = {m.v};

    bool contains_only_boundaries = m.v->is_boundary_vertex;
#ifdef MEMORY_DEBUG
    std::cout << "In expansion of " << print_v(m.v) << ", " << print_v(m.w) << ":";
    for (sptr<vertex_t> v : m.path) std::cout << " " << print_v(v);
    std::cout << std::endl;
#endif
    for (size_t i = 1; i < m.path.size(); i++) {
        sptr<vertex_t> v = m.path.at(i-1),
                        w = m.path.at(i);
        if (decoding_graph->share_hyperedge({v, w})) {
            curr.w = w;
            curr.path.push_back(w);
            contains_only_boundaries &= w->is_boundary_vertex;
        } else {
            sptr<hyperedge_t> e = get_flag_edge_for({v, w});
            error_chain_t ec = decoding_graph->get_error_chain(v, w, m.c1, m.c2, true);
            if (e != nullptr && (flag_ctr_map.at(e) == 2 || ec.path.size() > 3)) {
                // Finish off the current assignment and ignore the flag edge.
                push_back_assignment(assign_arr, curr);
#ifdef MEMORY_DEBUG
                std::cout << "\tFound flag edge [ " << print_v(v) << ", " << print_v(w) << " ], path:";
                for (sptr<vertex_t> x : ec.path) std::cout << " " << print_v(x);
                std::cout << std::endl;
                std::cout << "\t[F] pushing back assignment " << print_v(curr.v) << " <---> "
                    << print_v(curr.w) << std::endl;
#endif
                curr.v = w;
                curr.w = nullptr;
                curr.path = {w};
            } else {
                // Expand the path between v and w.
                if (ec.path.front() != v) std::reverse(ec.path.begin(), ec.path.end());
                for (size_t j = 1; j < ec.path.size(); j++) {
                    sptr<vertex_t> x = ec.path[j-1],
                                    y = ec.path[j];
                    curr.w = y;
                    curr.path.push_back(y);
                    contains_only_boundaries &= y->is_boundary_vertex;
                    if (y->is_boundary_vertex) {
                        if (!contains_only_boundaries) {
                            push_back_assignment(assign_arr, curr);
                        }
                        curr.v = y;
                        curr.w = nullptr;
                        curr.path = {y};
                        contains_only_boundaries = true;
                    }
                }
            }
        }
        if (w->is_boundary_vertex) {
            if (!contains_only_boundaries) {
                push_back_assignment(assign_arr, curr);
            }
            curr.v = w;
            curr.w = nullptr;
            curr.path = {w};
            contains_only_boundaries = true;
        }
    }
    if (!contains_only_boundaries) {
        push_back_assignment(assign_arr, curr);
    }
}

std::vector<component_t>
RestrictionDecoder::compute_connected_components(const std::vector<assign_t>& assignments) {
    // For each connected component, we know the following:
    //  (1) A boundary may be connected to multiple vertices (direct matches
    //      to the boundary, or matches that go through a boundary).
    //  (2) Non-boundary vertices are only connected to R-1 vertices (at most),
    //      where R is the number of restricted lattices.
    struct e_t : base::edge_t {
        assign_t m;
    };
    auto cgr = std::make_unique<Graph<vertex_t, e_t>>();
    // Add all assignments to the graph.
    for (const assign_t& m : assignments) {
        if (m.v == nullptr) continue;
        if (!cgr->contains(m.v)) cgr->add_vertex(m.v);
        if (!cgr->contains(m.w)) cgr->add_vertex(m.w);
        if (cgr->contains(m.v, m.w)) continue;
        sptr<e_t> e = cgr->make_and_add_edge(m.v, m.w);
        e->m = m;
    }
    // If the red boundary is not in the graph, then just exit as there will be
    // no connected components.
    sptr<vertex_t> vrb = decoding_graph->get_boundary_vertex(COLOR_RED);
    if (!cgr->contains(vrb)) return {};
    std::vector<component_t> components;

    typedef std::tuple<sptr<vertex_t>, std::vector<sptr<vertex_t>>, bool> c_entry_t;
    std::deque<c_entry_t> bfs;
    // Populate the data structures for adj(vrb).
    for (sptr<vertex_t> v : cgr->get_neighbors(vrb)) {
        bfs.emplace_back(v, std::vector<sptr<vertex_t>>{vrb, v}, true);
    }
    while (bfs.size()) {
        auto [v, path, is_entry_from_vrb] = bfs.front();
        bfs.pop_front();
        if (v->is_boundary_vertex) {
            // Then compute the connected component.
            std::vector<assign_t> assign_list;
            for (size_t i = 1; i < path.size(); i++) {
                sptr<vertex_t> x = path.at(i-1),
                                y = path.at(i);
                auto e = cgr->get_edge(x, y);
                assign_list.push_back(e->m);
            }
            // Make sure this is not identical to any other components.
            bool found_copy = false;
            std::vector<assign_t> rev_assign_list(assign_list.rbegin(), assign_list.rend());
            for (const auto& cc : components) {
                if (assign_list == cc.assignments || rev_assign_list == cc.assignments) {
                    found_copy = true;
                    break;
                }
            }
            if (found_copy) continue;
            int cc_color = get_complementary_colors_to(
                                {vrb->color, v->color}, decoding_graph->number_of_colors)[0];
            components.push_back({assign_list, cc_color});
#ifdef MEMORY_DEBUG
            std::cout << "Connected component:";
            for (sptr<vertex_t> x : path) std::cout << " " << print_v(x);
            std::cout << std::endl;
#endif
            continue;
        }
        // Move to the neighbors of v. If any of them are visited, then it implies we have a loop
        // to vrb.
        for (sptr<vertex_t> w : cgr->get_neighbors(v)) {
            if (w == vrb) {
                if (is_entry_from_vrb) continue;
            } else {
                if (std::find(path.begin(), path.end(), w) != path.end()) continue;
            }
            std::vector<sptr<vertex_t>> _path(path);
            _path.push_back(w);
            bfs.emplace_back(w, _path, false);
        }
    }
    return components;
}

std::set<vpair_t>
RestrictionDecoder::insert_error_chain_into(
        std::map<vpair_t, size_t>& incidence_map,
        const std::vector<sptr<vertex_t>>& path,
        int component_color,
        int c1,
        int c2)
{
    std::set<vpair_t> edge_set;
    for (size_t i = 1; i < path.size(); i++) {
        sptr<vertex_t> v = path.at(i-1),
                       w = path.at(i);
        if (v->is_boundary_vertex && w->is_boundary_vertex) { 
            continue;
        }
        // Flatten v and w.
        sptr<vertex_t> fv = v->get_base(),
                       fw = w->get_base();
        if (fv == fw) continue;
        // Make sure that their colors are not equal. Otherwise, get a common neighbor
        // of fv and fw which has a different color (call this fu). 
        // Add (fv, fu) and (fu, fw) instead.
        if (!decoding_graph->share_hyperedge({fv, fw})) {
            std::cerr << "not supposed to happen anymore..." << std::endl;
            exit(1);
        } else {
            if (fv->color != component_color && fw->color != component_color) {
                continue;
            }
            vpair_t e = make_vpair(fv, fw);
            incidence_map[e]++;
            edge_set.insert(e);
        }
    }
    return edge_set;
}

fp_t
RestrictionDecoder::lifting(
        stim::simd_bits_range_ref<SIMD_WIDTH> corr,
        const std::map<sptr<hyperedge_t>, sptr<hyperedge_t>>& best_rep_map,
        size_t tr) 
{
    fp_t out_log_pr = 0.0;

    const size_t MAX_TRIES = 31;
    const int color_off = tr/5;
    // Recursively call the lifting procedure until finished or tr is too large.
    //
    // Create incident vertex sets.
    std::set<sptr<vertex_t>> not_cc_incident, in_cc_incident;
    insert_incident_vertices(not_cc_incident, not_cc_map, color_plus_offset(COLOR_RED, color_off));
    for (auto& [es, c] : component_edge_sets) {
        for (auto it = es.begin(); it != es.end(); ) {
            if (!in_cc_map.count(*it)) it = es.erase(it);
            else                       it++;
        }
        insert_incident_vertices(in_cc_incident, es, color_plus_offset(c, color_off));
    }
    std::set<sptr<vertex_t>> all_incident(not_cc_incident);
    vtils::insert_range(all_incident, in_cc_incident);

#ifdef MEMORY_DEBUG
    std::cout << "Edges in CC:" << std::endl;
    for (auto& [e, cnt] : in_cc_map) {
        std::cout << "\t[ " << print_v(e.first) << " "
            << print_v(e.second) << " ], count = " << cnt << std::endl;
    }
    std::cout << "Edges not in CC:" << std::endl;
    for (auto& [e, cnt] : not_cc_map) {
        std::cout << "\t[ " << print_v(e.first) << " "
            << print_v(e.second) << " ], count = " << cnt << std::endl;
    }
#endif

    for (sptr<vertex_t> v : all_incident) {
        std::set<face_t> faces = get_faces(v, best_rep_map);
        const size_t nf = faces.size();
        const uint64_t enf = 1L << nf;

#ifdef MEMORY_DEBUG
        std::cout << "Faces of " << print_v(v) << ":" << std::endl;
        for (face_t fc : faces) {
            std::cout << "\t<";
            for (sptr<vertex_t> x : fc.vertices) std::cout << " " << print_v(x);
            std::cout << " >, frames =";
            for (uint64_t fr : fc.frames) std::cout << " " << fr;
            std::cout << std::endl;
        }
#endif
        // Track intersections with connected components and outside of
        // connected components.
        std::set<vpair_t> best_cc_boundary,
                          best_no_cc_boundary;
        stim::simd_bits<SIMD_WIDTH> best_cc_corr(corr.num_bits_padded()),
                                    best_no_cc_corr(corr.num_bits_padded());
        fp_t best_log_prob_cc = std::numeric_limits<fp_t>::lowest(),
             best_log_prob_no_cc = std::numeric_limits<fp_t>::lowest();
        for (uint64_t i = 0; i < enf; i++) {
            std::set<vpair_t> boundary;
            stim::simd_bits<SIMD_WIDTH> local_corr(corr.num_bits_padded());
            fp_t log_pr = 0.0;

            uint64_t ii = i;
            for (auto it = faces.begin(); it != faces.end() && ii; it++) {
                if (ii & 1) {
                    intersect_with_boundary(boundary, local_corr, log_pr, *it, v);
                }
                ii >>= 1;
            }
            // Now, we need to check how much boundary intersects with connected components
            // or anything not in the connected components. Note that this intersection is
            // in the neighborhood of v.
            size_t int_in_cc = locally_matches(in_cc_map, boundary, v),
                   int_not_cc = locally_matches(not_cc_map, boundary, v);
            update_best_boundary(
                    int_in_cc,
                    best_log_prob_cc,
                    best_cc_boundary,
                    best_cc_corr, 
                    log_pr,
                    boundary,
                    local_corr);
            update_best_boundary(
                    int_not_cc,
                    best_log_prob_no_cc,
                    best_no_cc_boundary,
                    best_no_cc_corr,
                    log_pr,
                    boundary,
                    local_corr);
        }
        if (best_cc_boundary.empty() && best_no_cc_boundary.empty()) continue;
        if (best_log_prob_cc > best_log_prob_no_cc) {
            update_correction(
                in_cc_map, corr, out_log_pr, best_cc_corr, best_cc_boundary, best_log_prob_cc);
#ifdef MEMORY_DEBUG
                std::cout << "Matched " << print_v(v) << " to CCs with boundary:";
                for (vpair_t e : best_cc_boundary) {
                    std::cout << " (" << print_v(e.first) << ", " << print_v(e.second) << ")";
                }
                std::cout << std::endl;
#endif
        } else {
            update_correction(
                not_cc_map, corr, out_log_pr, best_no_cc_corr, best_no_cc_boundary, best_log_prob_no_cc);
#ifdef MEMORY_DEBUG
            std::cout << "Matched " << print_v(v) << " to bulk with boundary:";
            for (vpair_t e : best_no_cc_boundary) {
                std::cout << " (" << print_v(e.first) << ", " << print_v(e.second) << ")";
            }
            std::cout << std::endl;
#endif
        }
    }
    // Remove any widowed edges, as these can cause the decoder to loop infinitely.
    if (remove_widowed_edges(in_cc_map) || remove_widowed_edges(not_cc_map)) {
        out_log_pr -= 100;
    }
    if (in_cc_map.size() > 1 || not_cc_map.size() > 1) {
        if (tr < MAX_TRIES) {
            return out_log_pr + lifting(corr, best_rep_map, tr+1);
        }
    }
    return out_log_pr;
}

std::map<sptr<hyperedge_t>, sptr<hyperedge_t>>
RestrictionDecoder::flatten_edge_map(const std::map<sptr<hyperedge_t>, sptr<hyperedge_t>>& edge_map) {
    std::map<sptr<hyperedge_t>, sptr<hyperedge_t>> flat_edge_map;
    for (const auto& [x, y] : edge_map) {
        sptr<hyperedge_t> fx = decoding_graph->get_base_edge(x),
                          fy = decoding_graph->get_base_edge(y);
        flat_edge_map[fx] = fy;
    }
    return flat_edge_map;
}

std::set<face_t>
RestrictionDecoder::get_faces(
        sptr<vertex_t> v,
        const std::map<sptr<hyperedge_t>, sptr<hyperedge_t>>& best_rep_map) 
{
    std::set<face_t> faces;
    // Here, we will flatten the faces themselves. We require that
    // each face has different colors.
    for (sptr<hyperedge_t> e : decoding_graph->get_common_hyperedges({v})) {
        e = best_rep_map.count(e) ? best_rep_map.at(e) : e;
        face_t fc = make_face(e);
        if (fc.vertices.size()) {
            faces.insert(make_face(e));
        }
    }
    return faces;
}

face_t
make_face(sptr<hyperedge_t> e) {
    face_t fc;
    fc.frames = e->frames;
    fc.probability = e->probability;
    // Now flatten all the vertices.
    std::set<int> colors_in_face;
    for (size_t i = 0; i < e->get_order(); i++) {
        sptr<vertex_t> fx = e->get<vertex_t>(i)->get_base();
        if (std::find(fc.vertices.begin(), fc.vertices.end(), fx) != fc.vertices.end()) {
            fc.vertices.clear();
            return fc;
        }
        if (colors_in_face.count(fx->color)) {
            fc.vertices.clear();
            return fc;
        }
        fc.vertices.push_back(fx);
        colors_in_face.insert(fx->color);
    }
    std::sort(fc.vertices.begin(), fc.vertices.end());
    return fc;
}

void
intersect_with_boundary(
            std::set<vpair_t>& boundary,
            stim::simd_bits_range_ref<SIMD_WIDTH> corr,
            fp_t& log_pr,
            const face_t& fc,
            sptr<vertex_t> v)
{
    // Get edges of the hyperedge.
    for (size_t j = 0; j < fc.vertices.size(); j++) {
        auto x = fc.vertices.at(j);
        for (size_t k = j+1; k < fc.vertices.size(); k++) {
            auto y = fc.vertices.at(k);
            // Make sure that one of x or y is v.
            if (x != v && y != v) continue;
            vpair_t xy = make_vpair(x, y);
            boundary ^= xy;
        }
    }
    for (uint64_t fr : fc.frames) corr[fr] ^= 1;
    log_pr += log(fc.probability);
}

}   // qontra
