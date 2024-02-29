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
        fp_t& best_p,
        std::set<vpair_t>& best_boundary,
        stim::simd_bits_range_ref<SIMD_WIDTH> best_corr,
        fp_t& p,
        std::set<vpair_t>& boundary,
        stim::simd_bits_range_ref<SIMD_WIDTH> corr)
{
    if (intersection > 0 && p > best_p) {
        best_p = p;
        best_boundary = boundary;
        best_corr.clear();
        best_corr |= corr;
        return true;
    }
    return false;
}

void
remove_widowed_edges(std::map<vpair_t, size_t>& incidence_map) {
    std::map<sptr<vertex_t>, size_t> vertex_inc_map;
    for (auto& p : incidence_map) {
        const vpair_t& e = p.first;
        sptr<vertex_t> v1 = e.first,
                       v2 = e.second;
        vertex_inc_map[v1]++;
        vertex_inc_map[v2]++;
    }
    // Remove any pairs of vertices where both endpoints only have a single
    // incidence.
    for (auto it = incidence_map.begin(); it != incidence_map.end(); ) {
        vpair_t e = it->first;
        sptr<vertex_t> v1 = e.first,
                       v2 = e.second;
        if (vertex_inc_map.at(v1) == 1 && vertex_inc_map.at(v2) == 1) {
            it = incidence_map.erase(it);
        } else {
            it++;
        }
    }
}

Decoder::result_t
RestrictionDecoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    const size_t n_obs = circuit.count_observables();
    stim::simd_bits<SIMD_WIDTH> corr(n_obs);

    load_syndrome(syndrome);

    auto _detectors = detectors;
    auto _flags = flags;

    std::cout << "syndrome: D[";
    for (uint64_t d : _detectors) std::cout << " " << d;
    std::cout << " ], F[";
    for (uint64_t f : _flags) std::cout << " " << f;
    std::cout << " ]" << std::endl;

    if (_detectors.empty()) return ret_no_detectors();

    corr ^= get_base_corr();

    // Compute the MWPM for each restricted lattice.
    std::vector<assign_t> matchings;
    for (int c1 = 0; c1 < decoding_graph->number_of_colors; c1++) {
        for (int c2 = c1+1; c2 < decoding_graph->number_of_colors; c2++) {
            load_syndrome(syndrome, c1, c2, false);
            std::vector<Decoder::assign_t> _matchings = compute_matching(c1, c2, true);

            std::cout << "Matchings on L(" << c1 << ", " << c2 << "):" << std::endl;
            for (Decoder::assign_t x : _matchings) {
                matchings.push_back(cast_assign(x, c1, c2));

                std::cout << "\t" << std::get<0>(x) << " <---> " << std::get<1>(x) << std::endl;
            }
        }
    }
    // Compute connected components.
    std::vector<component_t> components = compute_connected_components(matchings);
    // Identify all edges inside a connected component (and which component it is), and
    // do the same for all edges outside a connected component. Furthermore, maintain a
    // counter for the number of times an edge appears in a connected component.
    std::map<vpair_t, size_t> in_cc_map, not_cc_map;
    std::vector<std::pair<std::set<vpair_t>, int>> component_edge_sets;
    // Track all triggered flag edges (need these to update the correction)
    std::set<sptr<hyperedge_t>> triggered_flag_edges;
    // Track all matches that are in a component.
    std::set<assign_t> in_cc_assignments;
    for (const component_t& cc : components) {
        int color = cc.color;
        std::set<vpair_t> edge_set;
        for (const assign_t& m : cc.assignments) {
            in_cc_assignments.insert(m);
            sptr<vertex_t> v = decoding_graph->get_vertex(std::get<0>(m)),
                           w = decoding_graph->get_vertex(std::get<1>(m));
            std::set<vpair_t> tmp = insert_error_chain_into(
                                        in_cc_map,
                                        v,
                                        w,
                                        color,
                                        std::get<2>(m),
                                        std::get<3>(m),
                                        false,
                                        triggered_flag_edges);
            vtils::insert_range(edge_set, tmp);
        }
        component_edge_sets.emplace_back(edge_set, color);
    }
    // Do not_cc now.
    for (const assign_t& m : matchings) {
        if (in_cc_assignments.count(m)) continue;
        sptr<vertex_t> v = decoding_graph->get_vertex(std::get<0>(m)),
                       w = decoding_graph->get_vertex(std::get<1>(m));
        if (std::get<2>(m) != COLOR_RED && std::get<3>(m) != COLOR_RED) continue;
        insert_error_chain_into(
                not_cc_map,
                v,
                w,
                COLOR_RED,
                std::get<2>(m),
                std::get<3>(m),
                false,
                triggered_flag_edges);
    }

    for (sptr<hyperedge_t> e : triggered_flag_edges) {
        std::cout << "Applying correction for edge D[";
        for (sptr<vertex_t> v : e->get<vertex_t>()) std::cout << " " << print_v(v);
        std::cout << ", frames =";
        for (uint64_t fr : e->frames) std::cout << " " << fr;
        std::cout << std::endl;

        for (uint64_t fr : e->frames) corr[fr] ^= 1;
    }

    if (in_cc_map.empty() && not_cc_map.empty()) {
        return { 0, corr };
    }
    // Note that now, all edges have been flattened. We should also flatten the flags
    // for when we get the incident faces.
    for (uint64_t& f : flags) {
        f = circuit.detector_base_map.at(f);
    }
    decoding_graph->activate_flags(flags);

    // Here, we will iterate multiple times and try to match faces as much as possible.
    size_t tries = 0;
r_compute_correction:
    std::set<sptr<vertex_t>> not_cc_incident,
                             in_cc_incident;
    insert_incident_vertices(not_cc_incident, not_cc_map, (COLOR_RED + tries/5) % 3);
    for (auto& es : component_edge_sets) {
        // Remove any deleted edges from the edge set.
        for (auto it = es.first.begin(); it != es.first.end(); ) {
            if (!in_cc_map.count(*it))  it = es.first.erase(it);
            else                        it++;
        }
        insert_incident_vertices(in_cc_incident, es.first, (es.second + tries/5) % 3);
    }
    std::set<sptr<vertex_t>> all_incident(not_cc_incident);
    vtils::insert_range(all_incident, in_cc_incident);

    for (sptr<vertex_t> v : all_incident) {
        std::set<sptr<hyperedge_t>> faces = get_faces(v);
        const size_t nf = faces.size();
        if (nf > 20) {
            std::cerr << "[ RestrictionDecoder ] found vertex " << print_v(v)
                << " with too many faces: " << nf << std::endl;
            for (sptr<hyperedge_t> e : faces) {
                std::cerr << "\t<";
                for (size_t i = 0; i < e->get_order(); i++) {
                    std::cerr << " " << print_v(e->get<vertex_t>(i));
                }
                std::cerr << " >" << std::endl;
            }
            exit(1);
        }
        const uint64_t enf = 1L << nf;
        // Track intersections with connected components and outside of
        // connected components.
        std::set<vpair_t> best_cc_boundary,
                          best_no_cc_boundary;
        stim::simd_bits<SIMD_WIDTH> best_cc_corr(n_obs),
                                    best_no_cc_corr(n_obs);
        fp_t best_prob_cc = 0.0,
             best_prob_no_cc = 0.0;
        for (uint64_t i = 0; i < enf; i++) {
            std::set<vpair_t> boundary;
            stim::simd_bits<SIMD_WIDTH> local_corr(n_obs);
            fp_t pr = 1.0;

            uint64_t ii = i;
            for (auto it = faces.begin(); it != faces.end() && ii; it++) {
                if (ii & 1) {
                    intersect_with_boundary(boundary, local_corr, pr, *it, v);
                }
                ii >>= 1;
            }
            // Now, we need to check how much boundary intersects with connected components
            // or anything not in the connected components. Note that this intersection is
            // in the neighborhood of v.
            size_t int_in_cc = locally_matches(in_cc_map, boundary, v),
                   int_not_cc = locally_matches(not_cc_map, boundary, v);
            update_best_boundary(
                int_in_cc, best_prob_cc, best_cc_boundary, best_cc_corr, pr, boundary, local_corr);
            update_best_boundary(
                int_not_cc, best_prob_no_cc, best_no_cc_boundary, best_no_cc_corr, pr, boundary, local_corr);
        }
        if (best_prob_cc > best_prob_no_cc) {
            corr ^= best_cc_corr;
            for (const vpair_t& e : best_cc_boundary) {
                if ((--in_cc_map[e]) == 0) in_cc_map.erase(e);
            }
        } else {
            corr ^= best_no_cc_corr;
            for (const vpair_t& e : best_no_cc_boundary) {
                if ((--not_cc_map[e]) == 0) not_cc_map.erase(e);
            }
        }
    }
    // Remove any widowed edges, as these can cause the decoder to loop infinitely.
    remove_widowed_edges(in_cc_map);
    remove_widowed_edges(not_cc_map);
    if (in_cc_map.size() > 1 || not_cc_map.size() > 1) {
        if (tries < 31) {
            tries++;
            goto r_compute_correction;
        } else {
            // Throw an error here.
            std::cerr << "[ RestrictionDecoder ] Failed to compute correction for syndrome: D[";
            for (uint64_t d : _detectors) std::cerr << " " << d;
            std::cerr << " ] and F[";
            for (uint64_t f : _flags) std::cerr << " " << f;
            std::cerr << " ]" << std::endl
                    << "\tEdges remaining in CC:";
            for (auto& p : in_cc_map) {
                vpair_t e = p.first;
                std::cerr << " " << print_v(e.first) << "," << print_v(e.second);
            }
            std::cerr << std::endl << "\tEdges remaining out of CC:";
            for (auto& p : not_cc_map) {
                vpair_t e = p.first;
                std::cerr << " " << print_v(e.first) << "," << print_v(e.second);
            }
            std::cerr << std::endl;
            std::cerr << "Faces for incident vertices:" << std::endl;
            for (sptr<vertex_t> v : all_incident) {
                std::set<sptr<hyperedge_t>> faces = get_faces(v);
                const size_t nf = faces.size();
                std::cerr << "\t" << print_v(v) << " (nf = " << nf << "):" << std::endl;
                for (sptr<hyperedge_t> e : faces) {
                    std::cerr << "\t\t<";
                    for (size_t i = 0; i < e->get_order(); i++) {
                        std::cerr << " " << print_v(e->get<vertex_t>(i));
                    }
                    std::cerr << " >" << std::endl;
                }
            }
        }
    }
    // Otherwise, we are done.
    return { 0.0, corr };
}

std::vector<RestrictionDecoder::component_t>
RestrictionDecoder::compute_connected_components(
        const std::vector<RestrictionDecoder::assign_t>& assignments)
{
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
    for (const assign_t& m : assignments) {
        uint64_t dv = std::get<0>(m),
                 dw = std::get<1>(m);
        sptr<vertex_t> v = decoding_graph->get_vertex(dv),
                       w = decoding_graph->get_vertex(dw);
        if (!cgr->contains(v)) cgr->add_vertex(v);
        if (!cgr->contains(w)) cgr->add_vertex(w);
        if (cgr->contains(v, w)) continue;
        sptr<e_t> e = cgr->make_and_add_edge(v, w);
        e->c1 = std::get<2>(m);
        e->c2 = std::get<3>(m);
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
        std::vector<assign_t> assign_list{ std::make_tuple(v->id, vrb->id, e->c1, e->c2) };
        
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
        bool force_unflagged,
        std::set<sptr<hyperedge_t>>& triggered_flag_edges)
{
    error_chain_t ec = decoding_graph->get(c1, c2, src, dst, force_unflagged);

    std::cout << "error chain btwn " << print_v(src) << " and " << print_v(dst) << ":";
    for (sptr<vertex_t> x : ec.path) {
        std::cout << " " << print_v(x);
    }
    std::cout << std::endl;

    for (size_t i = 1; i < ec.path.size(); i++) {
        sptr<vertex_t> v = ec.path.at(i-1),
                       w = ec.path.at(i);
        if (v->is_boundary_vertex && w->is_boundary_vertex) { 
            continue;
        }
        // Update the correction as this maybe a flag edge.
        sptr<hyperedge_t> e = get_flag_edge_for({v, w});
        if (e != nullptr) {
            triggered_flag_edges.insert(e);
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
            sptr<vertex_t> fu = nullptr;
            for (sptr<vertex_t> x : decoding_graph->get_common_neighbors({fv, fw})) {
                sptr<vertex_t> fx = x->get_base();
                if (fx->color != c1 && fx->color != c2) continue;
                if (fx == fv || fx == fw) continue;
                if (fv->color != component_color 
                        && fw->color != component_color
                        && fx->color != component_color)
                {
                    continue;
                }
                fu = fx;
                break;
            }
            if (fu == nullptr) {
                // When this happens, just find a path in the non-flagged graph to use.
                insert_error_chain_into(
                        incidence_map, fv, fw, component_color, c1, c2, true, triggered_flag_edges);
            } else {
                vpair_t e1 = make_vpair(fv, fu),
                        e2 = make_vpair(fu, fw);
                incidence_map[e1]++;
                incidence_map[e2]++;
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
    for (auto& p : incidence_map) edge_set.insert(p.first);
    return edge_set;
}

std::set<sptr<hyperedge_t>>
RestrictionDecoder::get_faces(sptr<vertex_t> v) {
    std::set<sptr<hyperedge_t>> faces;
    // Here, we will flatten the faces themselves. We require that
    // each face has different colors.
    for (sptr<hyperedge_t> e : decoding_graph->get_common_hyperedges({v})) {
        std::vector<sptr<vertex_t>> flat_vlist;
        std::set<int> colors_in_face;
        bool do_not_add = false;
        for (size_t i = 0; i < e->get_order(); i++) {
            sptr<vertex_t> x = e->get<vertex_t>(i);
            sptr<vertex_t> fx = x->get_base();
            if (std::find(flat_vlist.begin(), flat_vlist.end(), fx) != flat_vlist.end()) {
                continue;
            }
            if (colors_in_face.count(fx->color)) {
                do_not_add = true;
                break;
            }
            flat_vlist.push_back(fx);
            colors_in_face.insert(fx->color);
        }
        if (do_not_add) continue;
        sptr<hyperedge_t> fe = decoding_graph->get_edge(flat_vlist);
        if (fe != nullptr) {
            fe = decoding_graph->get_best_edge_from_class_of(fe);
            faces.insert(fe);
        }
    }
    return faces;
}

}   // qontra
