/* 
 *  author: Suhas Vittal
 *  date:   17 February 2024
 * */

#include "qontra/decoder/restriction.h"

#include <vtils/set_algebra.h>
#include <vtils/utility.h>

#include <initializer_list>

const int COLOR_RED = 0;

namespace qontra {

using namespace graph;
using namespace decoding;

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

void
remove_widowed_edges(std::map<vpair_t, size_t>& incidence_map) {
    std::map<sptr<vertex_t>, size_t> vertex_inc_map;
    for (const auto& [e, cnt] : incidence_map) {
        const auto& [v1, v2] = e;
        vertex_inc_map[v1]++;
        vertex_inc_map[v2]++;
    }
    // Remove any pairs of vertices where both endpoints only have a single
    // incidence.
    for (auto it = incidence_map.begin(); it != incidence_map.end(); ) {
        const auto& [v1, v2] = it->first;
        if (vertex_inc_map.at(v1) == 1 && vertex_inc_map.at(v2) == 1) {
            it = incidence_map.erase(it);
        } else {
            it++;
        }
    }
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
    auto _detectors = detectors;
    auto _flags = flags;

    /*
    std::cout << "syndrome: D[";
    for (uint64_t d : _detectors) std::cout << " " << d;
    std::cout << " ], F[";
    for (uint64_t f : _flags) std::cout << " " << f;
    std::cout << " ]" << std::endl;
    */

    if (_detectors.empty()) return ret_no_detectors();

    corr ^= get_base_corr();

    // Compute the MWPM for each restricted lattice.
#ifdef DECODER_PERF
    timer.clk_start();
    fp_t t;
#endif
    std::vector<c_assign_t> matchings;
    for (int c1 = 0; c1 < decoding_graph->number_of_colors; c1++) {
        for (int c2 = c1+1; c2 < decoding_graph->number_of_colors; c2++) {
            load_syndrome(syndrome, c1, c2, false);
            std::vector<Decoder::assign_t> _matchings = compute_matching(c1, c2, true);

//          std::cout << "Matchings on L(" << c1 << ", " << c2 << "):" << std::endl;
            for (Decoder::assign_t x : _matchings) {
                matchings.push_back(cast_assign(x, c1, c2));

//              std::cout << "\t" << std::get<0>(x) << " <---> " << std::get<1>(x) << std::endl;
            }
        }
    }
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
    std::set<c_assign_t> in_cc_assignments;
    for (const component_t& cc : components) {
        int color = cc.color;
        std::set<vpair_t> edge_set;
        for (const c_assign_t& m : cc.assignments) {
            const auto& [d1, d2, c1, c2] = m;

            in_cc_assignments.insert(m);
            sptr<vertex_t> v = decoding_graph->get_vertex(d1),
                           w = decoding_graph->get_vertex(d2);
            std::set<vpair_t> tmp = insert_error_chain_into(in_cc_map, v, w, color, c1, c2, false);
            vtils::insert_range(edge_set, tmp);
        }
        component_edge_sets.emplace_back(edge_set, color);
    }
    // Do not_cc now.
    for (const c_assign_t& m : matchings) {
        if (in_cc_assignments.count(m)) continue;
        const auto& [d1, d2, c1, c2] = m;
        sptr<vertex_t> v = decoding_graph->get_vertex(d1),
                       w = decoding_graph->get_vertex(d2);
        if (c1 != COLOR_RED && c2 != COLOR_RED) continue;
        insert_error_chain_into(not_cc_map, v, w, COLOR_RED, c1, c2, false);
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
    std::cout << "corr1 = ";
    for (size_t i = 0; i < n_obs; i++) std::cout << corr1[i]+0;
    std::cout << std::endl;
    // If there are no triggered flag edges. Return here.
    if (triggered_flag_edges.empty()) {
        return { 0.0, corr1 };
    }
    std::cout << "Triggered flag edges:" << std::endl;
    for (auto& [e, path] : triggered_flag_edges) {
        std::cout << "\t[";
        for (sptr<vertex_t> v : e->get<vertex_t>()) std::cout << " " << print_v(v->get_base());
        std::cout << " ], frames =";
        for (uint64_t fr : e->frames) std::cout << " " << fr;
        std::cout << ", path = [";
        for (sptr<vertex_t> v : path) std::cout << " " << print_v(v->get_base());
        std::cout << " ]" << std::endl;
    }
    // Otherwise, perform the lifting procedure again, this time removing all edges corresponding
    // to triggered flag edges. Also update corr for each removed flag edge.
    in_cc_map = std::move(_in_cc_map);
    not_cc_map = std::move(_not_cc_map);
    component_edge_sets = std::move(_c_edge_sets);
    fp_t log_p2 = 0.0;

    std::set<sptr<hyperedge_t>> visited_flag_edges;
    for (const auto& [he, vlist] : triggered_flag_edges) {
        for (size_t i = 1; i < vlist.size(); i++) {
            sptr<vertex_t> v = vlist.at(i-1)->get_base(),
                           w = vlist.at(i)->get_base();
            vpair_t e = make_vpair(v, w);
            erase_from_incidence_map(e, in_cc_map);
            erase_from_incidence_map(e, not_cc_map);
        }
        // Check if there is a similar hyperedge.
        bool found_similar = false;
        for (sptr<hyperedge_t> _he : visited_flag_edges) {
            if (he->endpoints == _he->endpoints && he->frames == _he->frames) {
                found_similar = true;
                break;
            }
        }
        // Apply the edge's frame changes.
        if (!found_similar) {
            for (uint64_t fr : he->frames) corr2[fr] ^= 1;
            visited_flag_edges.insert(he);
        }
    }
    log_p2 += lifting(corr2, best_rep_map);
    std::cout << "corr1 = ";
    for (size_t i = 0; i < n_obs; i++) std::cout << corr1[i]+0;
    std::cout << std::endl;

    std::cout << "corr2 = ";
    for (size_t i = 0; i < n_obs; i++) std::cout << corr2[i]+0;
    std::cout << std::endl;
//  corr = std::move(log_p1 > log_p2 ? corr1 : corr2);
    return { 0.0, corr2 };
}

std::vector<component_t>
RestrictionDecoder::compute_connected_components(const std::vector<c_assign_t>& assignments) {
    // For each connected component, we know the following:
    //  (1) A boundary may be connected to multiple vertices (direct matches
    //      to the boundary, or matches that go through a boundary).
    //  (2) Non-boundary vertices are only connected to R-1 vertices (at most),
    //      where R is the number of restricted lattices.
    struct e_t : base::edge_t {
        int c1;
        int c2;
    };
    auto cgr = std::make_unique<Graph<vertex_t, e_t>>();
    // Add all assignments to the graph.
    for (const c_assign_t& m : assignments) {
        const auto& [dv, dw, c1, c2] = m;
        sptr<vertex_t> v = decoding_graph->get_vertex(dv),
                       w = decoding_graph->get_vertex(dw);
        if (!cgr->contains(v)) cgr->add_vertex(v);
        if (!cgr->contains(w)) cgr->add_vertex(w);
        if (cgr->contains(v, w)) continue;
        sptr<e_t> e = cgr->make_and_add_edge(v, w);
        e->c1 = c1;
        e->c2 = c2;
    }
    // If the red boundary is not in the graph, then just exit as there will be
    // no connected components.
    sptr<vertex_t> vrb = decoding_graph->get_boundary_vertex(COLOR_RED);
    if (!cgr->contains(vrb)) return {};
    std::vector<component_t> components;
    std::set<sptr<vertex_t>> skip_set;
    for (sptr<vertex_t> v : cgr->get_neighbors(vrb)) {
        if (skip_set.count(v)) continue;
        sptr<e_t> e = cgr->get_edge(v, vrb);
        std::vector<c_assign_t> assign_list{ std::make_tuple(v->id, vrb->id, e->c1, e->c2) };
        
        std::map<sptr<vertex_t>, sptr<vertex_t>> prev;
        prev[v] = vrb;
        sptr<vertex_t> curr = v;
        while (!curr->is_boundary_vertex) {
            // Compute next neighbor (should only be one).
            sptr<vertex_t> next = nullptr;
            for (sptr<vertex_t> w : cgr->get_neighbors(curr)) {
                if (w != prev[curr]) next = w;
            }
            if (next == nullptr) break;
            sptr<e_t> _e = cgr->get_edge(curr, next);
            uint64_t id1 = curr->id,
                     id2 = next->id;
            if (id1 > id2) {
                std::swap(id1, id2);
            }
            assign_list.emplace_back(id1, id2, _e->c1, _e->c2);
            prev[next] = curr;
            curr = next;
        }
        // If curr is a boundary, then we have discovered a connected component.
        if (curr->is_boundary_vertex) {
            // Note that if ended up at vrb anyways, we need to mark prev[curr] (prev[vrb])
            // to be skipped to avoid double counting.
            if (curr == vrb) skip_set.insert(prev[curr]);
            int cc_color = get_complementary_colors_to(
                                {vrb->color, curr->color}, decoding_graph->number_of_colors)[0];
            components.push_back({assign_list, cc_color});
        }
    }
    return components;
}

std::set<vpair_t>
RestrictionDecoder::insert_error_chain_into(
        std::map<vpair_t, size_t>& incidence_map,
        sptr<vertex_t> src,
        sptr<vertex_t> dst,
        int component_color,
        int c1,
        int c2,
        bool force_unflagged)
{
    error_chain_t ec = decoding_graph->get_error_chain(src, dst, c1, c2, force_unflagged);
    for (size_t i = 1; i < ec.path.size(); i++) {
        sptr<vertex_t> v = ec.path.at(i-1),
                       w = ec.path.at(i);
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
                // When this happens, just find a path in the non-flagged graph to use.
                insert_error_chain_into(incidence_map, v, w, component_color, c1, c2, true);
                // Update triggered flag edges if this is a flag edge.
                sptr<hyperedge_t> e = get_flag_edge_for({v, w});
                if (e != nullptr) {
                    error_chain_t ec = decoding_graph->get_error_chain(v, w, c1, c2, true);
                    triggered_flag_edges.emplace_back(e, ec.path);
                }
        } else {
            if (fv->color != component_color && fw->color != component_color) {
                continue;
            }
            vpair_t e = make_vpair(fv, fw);
            incidence_map[e]++;
        }
    }
    std::set<vpair_t> edge_set;
    for (auto& [e, cnt] : incidence_map) edge_set.insert(e);
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

    for (sptr<vertex_t> v : all_incident) {
        std::set<face_t> faces = get_faces(v, best_rep_map);
        const size_t nf = faces.size();
        const uint64_t enf = 1L << nf;
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
        if (best_log_prob_cc > best_log_prob_no_cc) {
            update_correction(
                in_cc_map, corr, out_log_pr, best_cc_corr, best_cc_boundary, best_log_prob_cc);
        } else {
            update_correction(
                not_cc_map, corr, out_log_pr, best_no_cc_corr, best_no_cc_boundary, best_log_prob_no_cc);
        }
    }
    // Remove any widowed edges, as these can cause the decoder to loop infinitely.
    remove_widowed_edges(in_cc_map);
    remove_widowed_edges(not_cc_map);
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
