/*
 *  author: Suhas Vittal
 *  date:   19 July 2023
 * */

#include "decoder/restriction.h"

namespace qontra {

//#define DEBUG

using namespace graph;
using namespace decoding;

template <class T> inline void
xor_entry_into(T x, std::set<T>& s) {
    if (s.count(x)) s.erase(x);
    else            s.insert(x);
}

inline std::string
print_v(colored_vertex_t* v) {
    std::string s;
    if (is_colored_boundary(v)) s += "B";
    else                s += std::to_string(v->id);
    s += "[" + v->color + "]";
    return s;
}

uint locally_matches(std::set<colored_edge_t*> s1, std::set<colored_edge_t*> s2, colored_vertex_t* incident) {
    std::set<colored_edge_t*> shared_edges;
    // We want to check that s2 is a subset of s1.
    for (auto e : s2) {
        if (e->src == incident || e->dst == incident) {
            if (!s1.count(e))   return 0;
            else                shared_edges.insert(e);
        }
    }
    return shared_edges.size();
}

Decoder::result_t
RestrictionDecoder::decode_error(stim::simd_bits_range_ref syndrome) {
    const uint det = circuit.count_detectors();
    const uint obs = circuit.count_observables();
#ifdef DEBUG
    std::cout << "===================================\n";
#endif

    std::vector<uint> detectors = get_nonzero_detectors(syndrome);
#ifdef DEBUG
    std::cout << "Detectors ( HW = " << detectors.size() << " ):";
    for (uint d : detectors) {
        auto v = c_decoding_graph.get_vertex(d);
        std::cout << " " << d << "[";
        if (circuit.flag_detection_events.count(d)) {
            std::cout << "F, " << v->color;
        } else {
            std::cout << (d % detectors_per_round) << ", " << v->color;
        }
        std::cout << "]";
    }
    std::cout << "\n";
    // Print flag stats.
    for (uint d : detectors) {
        if (!circuit.flag_detection_events.count(d)) continue;
        std::cout << d << "[F] hits:";
        auto v = c_decoding_graph.get_vertex(d);
        for (auto w : c_decoding_graph.get_neighbors(v)) {
            std::cout << " " << print_v(w);
            if (circuit.flag_detection_events.count(w->id)) std::cout << "[F]";
        }
        std::cout << "\n";
    }
#endif
    stim::simd_bits obs_bits(obs);
    for (uint i = 0; i < obs; i++) obs_bits[i] = syndrome[det+i];

    std::vector<match_t> matchings = blossom_subroutine(detectors);
    std::vector<cc_t> comps = compute_connected_components(matchings);
    // Now, divide the edges between those in the connected components
    // and those not in the components.
    typedef std::pair<std::set<colored_edge_t*>, std::string> edge_set_t;

    std::vector<edge_set_t> in_cc_array;

    std::set<colored_edge_t*> in_cc_set, not_cc_set;
    std::map<colored_edge_t*, uint> in_cc_count_map, not_cc_count_map;

    // Compute in cc.
    // 
    // We have to compute all edges in each cc (and out of all ccs later), and
    // we must do this up-to flattening.
    std::set<match_t> cc_match_set;
    for (cc_t& comp : comps) {
        std::vector<match_t> matches = std::get<2>(comp);
        std::set<colored_edge_t*> es;
        std::string cc_color = std::get<3>(comp);
#ifdef DEBUG
        std::cout << "connected component:\n";
#endif
        for (match_t& m : matches) {
            cc_match_set.insert(m);
            auto v = c_decoding_graph.get_vertex(std::get<0>(m));
            auto w = c_decoding_graph.get_vertex(std::get<1>(m));
#ifdef DEBUG
            std::cout << "\t" << print_v(v) << " <--> " << print_v(w) << ":";
#endif
            std::string r = std::get<2>(m);
            auto error_data = c_decoding_graph[r].get_error_chain_data(v, w);
            auto error_chain = error_data.error_chain;
            for (uint j = 1; j < error_chain.size(); j++) {
                colored_vertex_t* u1 = (colored_vertex_t*)error_chain[j-1];
                colored_vertex_t* u2 = (colored_vertex_t*)error_chain[j];
                if (is_colored_boundary(u1) && is_colored_boundary(u2)) continue;
                if (u1->color != cc_color && u2->color != cc_color) continue;
                // Flatten the vertices.
                colored_vertex_t* fu1 = flatten(u1);
                colored_vertex_t* fu2 = flatten(u2);
                auto e = c_decoding_graph.get_edge(fu1, fu2);
                if (e == nullptr) {
#ifdef DEBUG
                    std::cout << " (dne)" << print_v(fu1) << ", " << print_v(fu2);
#endif
                    continue;
                }
#ifdef DEBUG
                std::cout << " " << print_v(fu1) << ", " << print_v(fu2);
#endif
                xor_entry_into(e, es);
            }
#ifdef DEBUG
            std::cout << "\n";
#endif
        }
        for (auto e : es) {
            in_cc_set.insert(e);
            in_cc_count_map[e]++;
        }
        in_cc_array.push_back(std::make_pair(es, cc_color));
    }
    // Now compute out of cc, which contains anything not in cc.
#ifdef DEBUG
    std::cout << "not in any connected component:\n";
#endif
    for (match_t& m : matchings) {
        if (cc_match_set.count(m)) continue;
        colored_vertex_t* v = c_decoding_graph.get_vertex(std::get<0>(m));
        colored_vertex_t* w = c_decoding_graph.get_vertex(std::get<1>(m));

        std::string r = std::get<2>(m);
        if (r == "gb" || r == "bg") continue;
        auto error_data = c_decoding_graph[r].get_error_chain_data(v, w);
        auto error_chain = error_data.error_chain;
#ifdef DEBUG
            std::cout << "\t" << print_v(v) << " <--> " << print_v(w) << ":";
#endif
        for (uint j = 1; j < error_chain.size(); j++) {
            colored_vertex_t* u1 = (colored_vertex_t*)error_chain[j-1];
            colored_vertex_t* u2 = (colored_vertex_t*)error_chain[j];
            if (is_colored_boundary(u1) && is_colored_boundary(u2)) continue;
            // Flatten the vertices.
            colored_vertex_t* fu1 = flatten(u1);
            colored_vertex_t* fu2 = flatten(u2);
            auto e = c_decoding_graph.get_edge(fu1, fu2);
            if (e == nullptr) {
#ifdef DEBUG
                std::cout << " (dne)" << print_v(fu1) << ", " << print_v(fu2);
#endif
                continue;
            }
#ifdef DEBUG
            std::cout << " " << print_v(fu1) << ", " << print_v(fu2);
#endif
            not_cc_set.insert(e);
            not_cc_count_map[e]++;
        }
#ifdef DEBUG
        std::cout << "\n";
#endif
    }

    if (in_cc_set.empty() && not_cc_set.empty()) {
        stim::simd_bits corr(obs);
        corr.clear();
        return (Decoder::result_t) {
            0.0,
            corr,
            is_error(corr, syndrome)
        };
    }
    // Finally, compute the correction.
    stim::simd_bits corr(obs);
    corr.clear();
    uint tries = 0;
    // We will iteratively select faces to apply corrections on until
    // not_cc_incident and in_cc_incident are empty.
r_compute_correction:
    std::set<colored_vertex_t*> not_cc_incident = c_decoding_graph.get_all_incident_vertices(not_cc_set, "r");
    std::set<colored_vertex_t*> in_cc_incident;
    for (edge_set_t& es : in_cc_array) {
        // We need to filter out these edge sets -- we may have discarded certain edges.
        for (auto it = es.first.begin(); it != es.first.end(); ) {
            if (!in_cc_set.count(*it))  it = es.first.erase(it);
            else                        it++;
        }
        auto incident = c_decoding_graph.get_all_incident_vertices(es.first, es.second);
        for (auto x : incident) in_cc_incident.insert(x);
    }
    std::set<colored_vertex_t*> all_incident;
    for (auto v : not_cc_incident)  all_incident.insert(v);
    for (auto v : in_cc_incident)   all_incident.insert(v);

#ifdef DEBUG
    std::cout << "Edges in cc, out of cc: " << in_cc_set.size() << ", " << not_cc_set.size() << "\n";
    std::cout << "In CC:";
    for (auto e : in_cc_set) std::cout << " (" 
                                    << print_v((colored_vertex_t*)e->src)
                                    << ", "
                                    << print_v((colored_vertex_t*)e->dst)
                                    << ")";
    std::cout << "\nNot CC:";
    for (auto e : not_cc_set) std::cout << " (" 
                                    << print_v((colored_vertex_t*)e->src)
                                    << ", "
                                    << print_v((colored_vertex_t*)e->dst)
                                    << ")";
    std::cout << "\tincidents: " << in_cc_incident.size() << ", " << not_cc_incident.size() << "\n";
    std::cout << "\tIn CC:";
    for (auto v : in_cc_incident) std::cout << " " << print_v(v);
    std::cout << "\n\tNot CC:";
    for (auto v : not_cc_incident) std::cout << " " << print_v(v);
    std::cout << "\n";
#endif
    for (auto v : all_incident) {
        // Can only match one of not_cc or in_cc.
        std::set<face_t> incident_faces = get_incident_faces(v);
        uint64_t nf = incident_faces.size();
        std::vector<stim::simd_bits> face_corr_list;
        for (face_t fc : incident_faces) {
            face_corr_list.push_back(get_correction_for_face(fc));
        }
#ifdef DEBUG
        std::cout << "nf of " << print_v(v) << " = " << nf << ":\n";
        for (face_t f : incident_faces) {
            std::cout << "\t< " << print_v(std::get<0>(f))
                << ", " << print_v(std::get<1>(f))
                << ", " << print_v(std::get<2>(f)) 
                << ", " << get_correction_for_face(f)[0] << " >\n";
        }
#endif
        uint64_t enf = 1L << nf;

        uint64_t faces_cc = std::numeric_limits<uint64_t>::max();
        uint64_t faces_no_cc = std::numeric_limits<uint64_t>::max();
        uint64_t best_intersect_with_cc = 0;
        uint64_t best_intersect_with_no_cc = 0;
        // We need to track intersections on both sides.
        std::set<colored_edge_t*> best_cc_boundary, best_no_cc_boundary;
        stim::simd_bits best_cc_corr(obs);
        stim::simd_bits best_no_cc_corr(obs);
        for (uint64_t i = 0; i < enf; i++) {
            // Interpret the bits of i as the faces we will examine.
            std::set<colored_edge_t*> f_boundary;
            uint64_t j = i;
            auto it = incident_faces.begin();

            stim::simd_bits local_corr(obs);
            local_corr.clear();
            uint ii = 0;

            uint64_t faces = 0;
            while (j) {
                if (j & 1) {
                    face_t f = *it;
                    colored_vertex_t* v1 = std::get<0>(f);
                    colored_vertex_t* v2 = std::get<1>(f);
                    colored_vertex_t* v3 = std::get<2>(f);

                    colored_vertex_t* vertex_list[] = { v1, v2, v3 };
                    for (int r = 0; r < 3; r++) {
                        auto vx = vertex_list[r];
                        for (int s = r+1; s < 3; s++) {
                            auto vy = vertex_list[s];
                            auto exy = c_decoding_graph.get_edge(vx, vy);
                            if (vx == v || vy == v) {
                                if (exy == nullptr) {
                                    std::cout << print_v(vx) << ", " << print_v(vy)
                                        << " does not exist.\n";
                                }
                                xor_entry_into(exy, f_boundary);
                            }
                        }
                    }
                    faces++;
                    local_corr ^= face_corr_list[ii];
                }
                j >>= 1;
                it++;
                ii++;
            }
            // Check how much the edges on the boundary of the
            // faces intersect with in_cc or not_cc.
            uint int_with_cc = locally_matches(in_cc_set, f_boundary, v);
            uint int_with_no_cc = locally_matches(not_cc_set, f_boundary, v);
            if (int_with_cc == 0 && int_with_no_cc == 0) continue;

            if (int_with_cc > best_intersect_with_cc
                || (int_with_cc == best_intersect_with_cc && faces < faces_cc)) 
            {
                best_intersect_with_cc = int_with_cc;
                faces_cc = faces;
                best_cc_boundary = f_boundary;
                best_cc_corr = local_corr;
#ifdef DEBUG
                std::cout << "\tbest in cc now face set " << i << " (corr = " << local_corr[0] << ", int = " << best_intersect_with_cc << ")\n";
#endif
            } 
            if (int_with_no_cc > best_intersect_with_no_cc
                || (int_with_no_cc == best_intersect_with_no_cc && faces < faces_no_cc)) 
            {
                best_intersect_with_no_cc = int_with_no_cc;
                faces_no_cc = faces;
                best_no_cc_boundary = f_boundary;
                best_no_cc_corr = local_corr;
#ifdef DEBUG
                std::cout << "\tbest no cc now face set " << i << " (corr = " << local_corr[0] << ", int = " << best_intersect_with_no_cc << ")\n";
#endif
            }
        }
        std::set<colored_edge_t*> best_boundary;
        stim::simd_bits best_corr(obs);
        // Choose boundary with minimum size
        if (best_intersect_with_cc > best_intersect_with_no_cc
            || (best_intersect_with_cc == best_intersect_with_no_cc && v->color != "r")) 
        {
            best_boundary = best_cc_boundary;
            best_corr = best_cc_corr;
            for (auto e : best_boundary) {
                if ((--in_cc_count_map[e]) == 0) in_cc_set.erase(e);
            }
        } else {
            best_boundary = best_no_cc_boundary;
            best_corr = best_no_cc_corr;
            for (auto e : best_boundary) {
                if ((--not_cc_count_map[e]) == 0) not_cc_set.erase(e);
            }
        }
#ifdef DEBUG
        std::cout << "final face boundary has correction: " << best_corr[0] << ", true error: " << obs_bits[0] << "\n";
#endif
        // Commit the correction for the boundary and erase the edges.
        corr ^= best_corr;
    }
    // We should now filter out in_cc_set and not_cc_set for any widowed edges (edges that
    // cannot possibly form a face with any other edge).
    std::map<colored_vertex_t*, uint> in_cc_vertex_incidence_map, not_cc_vertex_incidence_map;
    for (auto e : in_cc_set) {
        in_cc_vertex_incidence_map[(colored_vertex_t*)e->src]++;
        in_cc_vertex_incidence_map[(colored_vertex_t*)e->dst]++;
    }
    for (auto e : not_cc_set) {
        not_cc_vertex_incidence_map[(colored_vertex_t*)e->src]++;
        not_cc_vertex_incidence_map[(colored_vertex_t*)e->dst]++;
    }
    // Now, we will just go remove any edge such that both endpoints have an incidence of one.
    for (auto it = in_cc_set.begin(); it != in_cc_set.end(); ) {
        colored_vertex_t* src = (colored_vertex_t*)(*it)->src;
        colored_vertex_t* dst = (colored_vertex_t*)(*it)->dst;
        if (in_cc_vertex_incidence_map[src] == 1 && in_cc_vertex_incidence_map[dst] == 1) {
            it = in_cc_set.erase(it);
        } else {
            it++;
        }
    }

    for (auto it = not_cc_set.begin(); it != not_cc_set.end(); ) {
        colored_vertex_t* src = (colored_vertex_t*)(*it)->src;
        colored_vertex_t* dst = (colored_vertex_t*)(*it)->dst;
        if (not_cc_vertex_incidence_map[src] == 1 && not_cc_vertex_incidence_map[dst] == 1) {
            it = not_cc_set.erase(it);
        } else {
            it++;
        }
    }

    if (in_cc_set.size() > 1 || not_cc_set.size() > 1) {
        if (tries < 100) {
            tries++;
            goto r_compute_correction;
        } else {
#ifdef DEBUG
            std::cout << "Failed to compute correction.\n";
            std::cout << "edges remaining in cc: " << in_cc_set.size() << "\n";
            for (auto e : in_cc_set) std::cout << "\t" << print_v((colored_vertex_t*)e->src) << ", " << print_v((colored_vertex_t*)e->dst) << "\n";
            std::cout << "edges remaining out of cc: " << not_cc_set.size() << "\n";
            for (auto e : not_cc_set) std::cout << "\t" << print_v((colored_vertex_t*)e->src) << ", " << print_v((colored_vertex_t*)e->dst) << "\n";
#endif
        }
    }
#ifdef DEBUG
    std::cout << "final correction: " << corr[0] << ", true error: " << obs_bits[0] << "\n";
    std::cout << "is error: " << is_error(corr, syndrome) << "\n";
#endif

    return (Decoder::result_t) {
        0.0,
        corr,
        is_error(corr, syndrome)
    };
}

std::vector<RestrictionDecoder::match_t>
RestrictionDecoder::blossom_subroutine(const std::vector<uint>& detectors) {
    std::set<vertex_t*> all_flags;
    for (uint df : circuit.flag_detection_events) all_flags.insert(c_decoding_graph.get_vertex(df));
    // Partition the detectors into syndromes for each restricted lattice and
    // compute the MWPM.
    std::array<std::vector<uint>, 3> restricted_syndromes_stabilizers;
    std::array<std::vector<uint>, 3> restricted_syndromes_flags;
    restricted_syndromes_stabilizers.fill(std::vector<uint>());
    restricted_syndromes_flags.fill(std::vector<uint>());
    uint n = 0;
    for (uint d : detectors) {
        auto v = c_decoding_graph.get_vertex(d);
        int i1, i2;
        if (v->color == "r") {
            i1 = 0; i2 = 1;
        } else if (v->color == "g") {
            i1 = 0; i2 = 2;
        } else {
            i1 = 1; i2 = 2;
        }
        
        if (circuit.flag_detection_events.count(d)) {
            restricted_syndromes_flags[i1].push_back(d);
            restricted_syndromes_flags[i2].push_back(d);
        } else {
            restricted_syndromes_stabilizers[i1].push_back(d);
            restricted_syndromes_stabilizers[i2].push_back(d);
        }
    }

    // Add boundaries wherever needed.
    for (auto& rs : restricted_syndromes_stabilizers) {
        if (rs.size() & 1) rs.push_back(BOUNDARY_INDEX);
    }

    // Perform MWPM. As this is a bit specialized, we just implement by hand.
    const std::string rcolors[] = { "rg", "rb", "gb" };
    std::vector<match_t> match_list;
    for (uint k = 0; k < 3; k++) {
        std::string rc = rcolors[k];
#ifdef DEBUG
        std::cout << "For restricted lattice " << k << ":\n";
#endif
        std::vector<uint> stab_events = restricted_syndromes_stabilizers[k],
                            flag_events = restricted_syndromes_flags[k];
        const uint n = stab_events.size();
        const uint m = (n*(n+1))/2;
        PerfectMatching pm(n, m);
        pm.options.verbose = false;

        // Create flagged decoding graph if necessary.
        const bool flags_are_active = !flag_events.empty();
        if (flags_are_active) {
            std::vector<vertex_t*> detector_list;
            std::set<vertex_t*> flag_set;
            for (auto d : stab_events) {
                if (d == BOUNDARY_INDEX) {
                    if (rc[0] == 'r' || rc[1] == 'r') {
                        detector_list.push_back(c_decoding_graph.get_vertex(RED_BOUNDARY_INDEX));
                    }
                    if (rc[0] == 'g' || rc[1] == 'g') {
                        detector_list.push_back(c_decoding_graph.get_vertex(GREEN_BOUNDARY_INDEX));
                    }
                    if (rc[0] == 'b' || rc[1] == 'b') {
                        detector_list.push_back(c_decoding_graph.get_vertex(BLUE_BOUNDARY_INDEX));
                    }
                } else {
                    detector_list.push_back(c_decoding_graph.get_vertex(d));
                }
            }
            for (auto f : flag_events) {
                flag_set.insert(c_decoding_graph.get_vertex(f));
            }
            c_decoding_graph[rc].setup_flagged_decoding_graph(detector_list, flag_set, all_flags);
        }
        
        std::map<colored_vertex_t*, colored_vertex_t*> node_to_pref_boundary;

#define __GET_ERROR_CHAIN(v, w) flags_are_active\
                                ? c_decoding_graph[rc].get_error_chain_data_from_flagged_graph((v), (w))\
                                : c_decoding_graph[rc].get_error_chain_data((v), (w))
        for (uint i = 0; i < n; i++) {
            uint di = stab_events[i];
            auto vi = c_decoding_graph.get_vertex(di);
            for (uint j = i+1; j < n; j++) {
                uint dj = stab_events[j];
                colored_vertex_t* vj;
                if (dj == BOUNDARY_INDEX) {
                    // We need to get the correct boundary, which is the one
                    // closer to vi.
                    uint64_t b1, b2;
                    if (rc[0] == 'r')   b1 = RED_BOUNDARY_INDEX;
                    else                b1 = GREEN_BOUNDARY_INDEX;

                    if (rc[1] == 'g')   b2 = GREEN_BOUNDARY_INDEX;
                    else                b2 = BLUE_BOUNDARY_INDEX;

                    auto vb1 = c_decoding_graph.get_vertex(b1);
                    auto vb2 = c_decoding_graph.get_vertex(b2);

                    auto b1_error_data = __GET_ERROR_CHAIN(vi, vb1);
                    auto b2_error_data = __GET_ERROR_CHAIN(vi, vb2);
                    if (b1_error_data.weight < b2_error_data.weight) {
                        vj = vb1;
                    } else {
                        vj = vb2;
                    }
                    node_to_pref_boundary[vi] = vj;
                } else {
                    vj = c_decoding_graph.get_vertex(dj);
                }
                auto error_data = __GET_ERROR_CHAIN(vi, vj);
                
                uint32_t edge_weight;
                if (error_data.weight > 1000.0) edge_weight = 1000000000;
                else                            edge_weight = (uint32_t) 1000.0*error_data.weight;
                pm.AddEdge(i, j, edge_weight);
            }
        }
        pm.Solve();
        for (uint i = 0; i < n; i++) {
            uint j = pm.GetMatch(i);
            if (i >= j) continue;
            uint di = stab_events[i],
                dj = stab_events[j];
            auto vi = c_decoding_graph.get_vertex(di);
            auto vj = dj == BOUNDARY_INDEX ? node_to_pref_boundary[vi] : c_decoding_graph.get_vertex(dj);
            // Split the match into two if it passes through the boundary.
            colored_vertex_t* bi;
            colored_vertex_t* bj;
            if (c_decoding_graph.are_matched_through_boundary(vi, vj, rcolors[k], &bi, &bj)) {
                match_t m1 = make_match(vi->id, bi->id, rcolors[k]);
                match_t m2 = make_match(vj->id, bj->id, rcolors[k]);
                match_list.push_back(m1);
                match_list.push_back(m2);
            } else {
                match_t m = make_match(vi->id, vj->id, rcolors[k]);
                match_list.push_back(m);
            }
#ifdef DEBUG
            auto error_data = __GET_ERROR_CHAIN(vi, vj);
            std::cout << "\t" << print_v(vi) << " <--> " << print_v(vj) << " (thru b = " 
                << error_data.error_chain_runs_through_boundary << ")\n";
#endif
        }
    }
    return match_list;
}

std::vector<RestrictionDecoder::cc_t>
RestrictionDecoder::compute_connected_components(const std::vector<RestrictionDecoder::match_t>& match_list) {
    // For each connected component, note the following:
    //  (a) The boundaries may be connected to multiple vertices.
    //  (b) Each non-boundary vertex is connected to at most vertices, but note that
    //  when considering paths from one boundary to another, a vertex is either
    //  a dead end (one unique connection) or has only one option (two
    //  connections).
    // 
    // Define a connectivity graph where A is connected to B if they are
    // matched. Furthermore, connect A to a boundary if A's matching to B
    // passes through a boundary.
    typedef colored_vertex_t    v_t;
    struct e_t : base::edge_t {
        std::string color;
    };
    typedef Graph<v_t, e_t>     cgr_t;  // Connectivity Graph type.

    cgr_t connectivity_graph;
    connectivity_graph.dealloc_on_delete = false;
    // Add all boundaries to the graph.
    uint64_t boundary_indices[] = {RED_BOUNDARY_INDEX, GREEN_BOUNDARY_INDEX, BLUE_BOUNDARY_INDEX};
    for (uint64_t b : boundary_indices) {
        auto vb = c_decoding_graph.get_vertex(b);
        connectivity_graph.add_vertex(vb);
    }
    for (match_t m : match_list) {
        auto d1 = std::get<0>(m);
        auto d2 = std::get<1>(m);
        auto r = std::get<2>(m);

        auto v1 = c_decoding_graph.get_vertex(d1);
        auto v2 = c_decoding_graph.get_vertex(d2);

        if (!connectivity_graph.contains(v1)) connectivity_graph.add_vertex(v1);
        if (!connectivity_graph.contains(v2)) connectivity_graph.add_vertex(v2);

        e_t* e = new e_t;
        e->src = (void*)v1;
        e->dst = (void*)v2;
        e->color = r;
        connectivity_graph.add_edge(e);
    }
    std::vector<cc_t> cc_list;
    // Compute connected components.
    auto rb = c_decoding_graph.get_vertex(RED_BOUNDARY_INDEX);  // Look from
    std::set<v_t*> skip_set;     // If we have a loop CC, then we need to skip certain searches.
    for (auto n : connectivity_graph.get_neighbors(rb)) {
        if (skip_set.count(n))  continue;
        e_t* rb_n_e = connectivity_graph.get_edge(rb, n);
        std::vector<match_t> match_list{make_match(RED_BOUNDARY_INDEX, n->id, rb_n_e->color)};

        std::map<v_t*, v_t*> prev;
        prev[n] = rb;
        auto curr = n;
        while (!is_colored_boundary(curr)) {
            // Get next neighbor.
            v_t* next = nullptr;
            for (auto w : connectivity_graph.get_neighbors(curr)) {
                if (w != prev[curr])    next = w;
            }
            if (next == nullptr)    break;
            prev[next] = curr;
            e_t* e = connectivity_graph.get_edge(curr, next);
            match_list.push_back(make_match(curr->id, next->id, e->color));
            curr = next;
        }
        // Only proceeded if we ended up at a boundary.
        if (is_colored_boundary(curr)) {
            uint64_t b1 = RED_BOUNDARY_INDEX;
            uint64_t b2 = curr->id;

            std::string cc_color;
            if (b2 == RED_BOUNDARY_INDEX || b2 == BLUE_BOUNDARY_INDEX) {
                cc_color = "g";
                if (b2 == RED_BOUNDARY_INDEX) {
                    skip_set.insert(prev[curr]);
                }
            } else if (b2 == GREEN_BOUNDARY_INDEX) {
                cc_color = "b";
            }
            cc_list.push_back(std::make_tuple(b1, b2, match_list, cc_color));
        }
    }
    // Delete edges in graph before returning.
    for (auto e : connectivity_graph.get_edges()) delete e;
    return cc_list;
}

colored_vertex_t*
RestrictionDecoder::flatten(colored_vertex_t* v) {
    if (is_colored_boundary(v)) return v;
    else {
        // In a circuit-level error model, flattening to the lowest level of detectors
        // (the first round) is not good as this ignores CNOT errors that occur in a prior
        // round. Thus, we instead flatten to the second round of detectors.
        uint64_t id = (v->id % detectors_per_round);
        colored_vertex_t* fv = c_decoding_graph.get_vertex(id);
        return fv;
    }
}

std::set<face_t>
RestrictionDecoder::get_incident_faces(colored_vertex_t* v) {
    // We will be accessing this often. Just memoize.
    static std::map<colored_vertex_t*, std::set<face_t>> memo;
    if (memo.count(v)) return memo[v];

    std::set<face_t> face_set;

    v = flatten(v);
    auto v_adj = c_decoding_graph.get_neighbors(v);
    for (auto w : v_adj) {
        w = flatten(w);
        if (v == w) continue;
        if (v->color == w->color) continue;
        if (!c_decoding_graph.contains(v, w)) continue;
        auto vw_adj = c_decoding_graph.get_common_neighbors(v, w);
        for (auto u : vw_adj) {
            u = flatten(u);
            if (u->color == v->color || u->color == w->color) continue;
            if (u == v || u == w) continue;
            if (!c_decoding_graph.contains(u, v) || !c_decoding_graph.contains(u, w)) continue;
            face_t fc = make_face(v, w, u);
            face_set.insert(fc);
        }
    }
    memo[v] = face_set;
    return face_set;
}

stim::simd_bits
RestrictionDecoder::get_correction_for_face(face_t fc) {
    stim::simd_bits corr(circuit.count_observables());
    corr.clear();
    for (uint fr : c_decoding_graph.get_face_frame_changes(fc)) corr[fr] ^= 1;
    return corr;
}

}   // qontra
