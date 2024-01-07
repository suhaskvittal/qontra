/*
 *  author: Suhas Vittal
 *  date:   19 July 2023
 * */

#include "decoder/restriction.h"

namespace qontra {

using namespace graph;
using namespace decoding;

#define DEBUG

template <class T> inline void
xor_entry_into(T x, std::set<T>& s) {
    if (s.count(x)) s.erase(x);
    else            s.insert(x);
}

template <typename PTR> std::string
print_v(PTR v) {
    auto _v = std::reinterpret_pointer_cast<colored_vertex_t>(v);

    std::string s;
    if (is_colored_boundary(_v))    s += "B";
    else                            s += std::to_string(_v->id);
    s += "[" + _v->color + "]";
    return s;
}

uint locally_matches(
        std::set<sptr<colored_edge_t>> s1, 
        std::set<sptr<colored_edge_t>> s2,
        sptr<colored_vertex_t> incident) 
{
    std::set<sptr<colored_edge_t>> shared_edges;
    // We want to check that s2 is a subset of s1.
    for (sptr<colored_edge_t> e : s2) {
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
    // Check if any flags are set.
    flags_are_active = false;
    for (uint d : detectors) {
        flags_are_active |= circuit.flag_detection_events.count(d);
    }
#ifdef DEBUG
    std::cout << "Detectors ( HW = " << detectors.size() << " ):";
    for (uint d : detectors) {
        sptr<colored_vertex_t> v = c_decoding_graph.get_vertex(d);
        std::cout << " " << d << "[";
        if (circuit.flag_detection_events.count(d)) {
            std::cout << "F, " << v->color;
        } else {
            std::cout << (d % detectors_per_round) << ", " << v->color;
        }
        std::cout << "]";
    }
    std::cout << "\n";
#endif
    stim::simd_bits obs_bits(obs);
    for (uint i = 0; i < obs; i++) obs_bits[i] = syndrome[det+i];

    std::vector<match_t> matchings = blossom_subroutine(detectors);
    std::vector<cc_t> comps = compute_connected_components(matchings);
    // Now, divide the edges between those in the connected components
    // and those not in the components.
    typedef std::pair<std::set<sptr<colored_edge_t>>, std::string> edge_set_t;

    std::vector<edge_set_t> in_cc_array;

    std::set<sptr<colored_edge_t>> in_cc_set, not_cc_set;
    std::map<sptr<colored_edge_t>, uint> in_cc_count_map, not_cc_count_map;

    // Compute in cc.
    // 
    // We have to compute all edges in each cc (and out of all ccs later), and
    // we must do this up-to flattening.
    std::set<match_t> cc_match_set;
    for (cc_t& comp : comps) {
        std::vector<match_t> matches = std::get<2>(comp);
        std::string cc_color = std::get<3>(comp);

        std::set<sptr<colored_edge_t>> es;
        for (match_t& m : matches) {
            cc_match_set.insert(m);
            sptr<colored_vertex_t> v = c_decoding_graph.get_vertex(std::get<0>(m)),
                                    w = c_decoding_graph.get_vertex(std::get<1>(m));
            std::string r = std::get<2>(m);
            insert_error_chain_into(es, in_cc_count_map, cc_color, v, w, r);
        }
        for (auto e : es) in_cc_set.insert(e);
        in_cc_array.push_back(std::make_pair(es, cc_color));
    }
    // Now compute out of cc, which contains anything not in cc.
#ifdef DEBUG
    std::cout << "not in any connected component:\n";
#endif
    for (match_t& m : matchings) {
        if (cc_match_set.count(m)) continue;
        sptr<colored_vertex_t> v = c_decoding_graph.get_vertex(std::get<0>(m));
        sptr<colored_vertex_t> w = c_decoding_graph.get_vertex(std::get<1>(m));

        std::string r = std::get<2>(m);
        if (r == "gb" || r == "bg") continue;
        insert_error_chain_into(not_cc_set, not_cc_count_map, "r", v, w, r);
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
    std::set<sptr<colored_vertex_t>> not_cc_incident = c_decoding_graph.get_all_incident_vertices(not_cc_set, "r");
    std::set<sptr<colored_vertex_t>> in_cc_incident;
    for (edge_set_t& es : in_cc_array) {
        // We need to filter out these edge sets -- we may have discarded certain edges.
        for (auto it = es.first.begin(); it != es.first.end(); ) {
            if (!in_cc_set.count(*it))  it = es.first.erase(it);
            else                        it++;
        }
        auto incident = c_decoding_graph.get_all_incident_vertices(es.first, es.second);
        for (auto x : incident) in_cc_incident.insert(x);
    }
    std::set<sptr<colored_vertex_t>> all_incident;
    for (auto v : not_cc_incident)  all_incident.insert(v);
    for (auto v : in_cc_incident)   all_incident.insert(v);

#ifdef DEBUG
    std::cout << "Edges in cc, out of cc: " << in_cc_set.size() << ", " << not_cc_set.size() << "\n";
    std::cout << "In CC:";
    for (auto e : in_cc_set) std::cout << " (" 
                                    << print_v(e->src)
                                    << ", "
                                    << print_v(e->dst)
                                    << ")";
    std::cout << "\nNot CC:";
    for (auto e : not_cc_set) std::cout << " (" 
                                    << print_v(e->src)
                                    << ", "
                                    << print_v(e->dst)
                                    << ")";
    std::cout << "\tincidents: " << in_cc_incident.size() << ", " << not_cc_incident.size() << "\n";
    std::cout << "\tIn CC:";
    for (auto v : in_cc_incident) std::cout << " " << print_v(v);
    std::cout << "\n\tNot CC:";
    for (auto v : not_cc_incident) std::cout << " " << print_v(v);
    std::cout << "\n";
#endif
    bool no_progress = true;
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

        // We need to track intersections on both sides.
        std::set<sptr<colored_edge_t>> best_cc_boundary, best_no_cc_boundary;
        stim::simd_bits best_cc_corr(obs);
        stim::simd_bits best_no_cc_corr(obs);
        fp_t best_prob_cc = 0.0;
        fp_t best_prob_no_cc = 0.0;
        for (uint64_t i = 0; i < enf; i++) {
            // Interpret the bits of i as the faces we will examine.
            std::set<sptr<colored_edge_t>> f_boundary;
            uint64_t j = i;
            auto it = incident_faces.begin();

            stim::simd_bits local_corr(obs);
            local_corr.clear();
            uint ii = 0;

            uint64_t faces = 0;
            fp_t pr = 1.0;
            while (j) {
                if (j & 1) {
                    face_t f = *it;
                    sptr<colored_vertex_t> v1 = std::get<0>(f);
                    sptr<colored_vertex_t> v2 = std::get<1>(f);
                    sptr<colored_vertex_t> v3 = std::get<2>(f);

                    sptr<colored_vertex_t> vertex_list[] = { v1, v2, v3 };
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
                    pr *= c_decoding_graph.get_face_probability(f);
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

            if (int_with_cc > 0 && pr > best_prob_cc) {
                best_prob_cc = pr;
                best_cc_boundary = f_boundary;
                best_cc_corr = local_corr;
#ifdef DEBUG
                std::cout << "\tbest in cc now face set " << i << " (corr = " << local_corr[0] << ", int = " << int_with_cc << ", pr = " << pr << ")\n";
#endif
            } 
            if (int_with_no_cc > 0 && pr > best_prob_no_cc) {
                best_prob_no_cc = pr;
                best_no_cc_boundary = f_boundary;
                best_no_cc_corr = local_corr;
#ifdef DEBUG
                std::cout << "\tbest no cc now face set " << i << " (corr = " << local_corr[0] << ", int = " << int_with_no_cc << ", pr = " << pr << ")\n";
#endif
            }
        }
        std::set<sptr<colored_edge_t>> best_boundary;
        stim::simd_bits best_corr(obs);
        // Choose boundary with minimum size
        if (best_prob_cc > best_prob_no_cc) {
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
        no_progress &= best_boundary.empty();
    }
    // We should now filter out in_cc_set and not_cc_set for any widowed edges (edges that
    // cannot possibly form a face with any other edge).
    std::map<sptr<colored_vertex_t>, uint> in_cc_vertex_incidence_map, not_cc_vertex_incidence_map;
    for (auto e : in_cc_set) {
        auto src = std::reinterpret_pointer_cast<colored_vertex_t>(e->src);
        auto dst = std::reinterpret_pointer_cast<colored_vertex_t>(e->dst);
        in_cc_vertex_incidence_map[src]++;
        in_cc_vertex_incidence_map[dst]++;
    }
    for (auto e : not_cc_set) {
        auto src = std::reinterpret_pointer_cast<colored_vertex_t>(e->src);
        auto dst = std::reinterpret_pointer_cast<colored_vertex_t>(e->dst);
        not_cc_vertex_incidence_map[src]++;
        not_cc_vertex_incidence_map[dst]++;
    }
    // Now, we will just go remove any edge such that both endpoints have an incidence of one.
    for (auto it = in_cc_set.begin(); it != in_cc_set.end(); ) {
        sptr<colored_vertex_t> src = std::reinterpret_pointer_cast<colored_vertex_t>((*it)->src);
        sptr<colored_vertex_t> dst = std::reinterpret_pointer_cast<colored_vertex_t>((*it)->dst);
        if (in_cc_vertex_incidence_map[src] == 1 && in_cc_vertex_incidence_map[dst] == 1) {
            it = in_cc_set.erase(it);
        } else {
            it++;
        }
    }

    for (auto it = not_cc_set.begin(); it != not_cc_set.end(); ) {
        sptr<colored_vertex_t> src = std::reinterpret_pointer_cast<colored_vertex_t>((*it)->src);
        sptr<colored_vertex_t> dst = std::reinterpret_pointer_cast<colored_vertex_t>((*it)->dst);
        if (not_cc_vertex_incidence_map[src] == 1 && not_cc_vertex_incidence_map[dst] == 1) {
            it = not_cc_set.erase(it);
        } else {
            it++;
        }
    }

    if (in_cc_set.size() > 1 || not_cc_set.size() > 1) {
        if (tries < 100) {
            tries++;
            if (no_progress) {
                if (tries < 10) {
                    // Change all boundaries in case of issues. Maybe wrong boundary.
                    switch_out_boundaries(in_cc_set, in_cc_count_map);
                    switch_out_boundaries(not_cc_set, not_cc_count_map);
                } else {
                    // It could also be some weird boundary CX edge.
                    remap_boundary_edges(in_cc_set, in_cc_count_map);
                    remap_boundary_edges(not_cc_set, not_cc_count_map);
                }
            }
            goto r_compute_correction;
        } else {
#ifdef DEBUG
            std::cout << "Failed to compute correction.\n";
            std::cout << "edges remaining in cc: " << in_cc_set.size() << "\n";
            for (auto e : in_cc_set) std::cout << "\t" << print_v(e->src) << ", " << print_v(e->dst) << "\n";
            std::cout << "edges remaining out of cc: " << not_cc_set.size() << "\n";
            for (auto e : not_cc_set) std::cout << "\t" << print_v(e->src) << ", " << print_v(e->dst) << "\n";
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
    std::set<sptr<vertex_t>> all_flags;
    for (uint df : circuit.flag_detection_events) all_flags.insert(c_decoding_graph.get_vertex(df));
    // Partition the detectors into syndromes for each restricted lattice and
    // compute the MWPM.
    std::array<std::vector<uint>, 3> restricted_syndromes_stabilizers;
    std::vector<uint> flag_events;

    restricted_syndromes_stabilizers.fill(std::vector<uint>());
    for (uint d : detectors) {
        if (circuit.flag_detection_events.count(d)) {
            flag_events.push_back(d);
        } else {
            auto v = c_decoding_graph.get_vertex(d);
            int i1, i2;
            if (v->color == "r") {
                i1 = 0; i2 = 1;
            } else if (v->color == "g") {
                i1 = 0; i2 = 2;
            } else {
                i1 = 1; i2 = 2;
            }
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
        std::vector<uint> stab_events = restricted_syndromes_stabilizers[k];
        const uint n = stab_events.size();
        const uint m = (n*(n+1))/2;
        PerfectMatching pm(n, m);
        pm.options.verbose = false;

        // Create flagged decoding graph if necessary.
        if (flags_are_active) {
            std::vector<sptr<vertex_t>> detector_list;
            std::vector<DecodingGraph::flag_edge_t> flag_edge_list;
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
                auto flag_edge = circuit.flag_edge_table[f];
                uint src_id = std::get<0>(flag_edge),
                     thru_id = std::get<1>(flag_edge),
                     dst_id = std::get<2>(flag_edge);
                auto vthru = c_decoding_graph.get_vertex(thru_id);
                // Need to handle the case where one of src or dst are the boundary.
                if (src_id == -1) std::swap(src_id, dst_id);
                auto vsrc = c_decoding_graph.get_vertex(src_id);

                // If the source is not even in the restricted lattice, continue.
                if (!c_decoding_graph[rc].contains(vsrc)) continue;

                sptr<colored_vertex_t> vdst;
                if (dst_id == BOUNDARY_INDEX) {
                    if (vsrc->color == "r") {
                        vdst = c_decoding_graph.get_vertex(RED_BOUNDARY_INDEX);
                    } else if (vsrc->color == "g") {
                        vdst = c_decoding_graph.get_vertex(GREEN_BOUNDARY_INDEX);
                    } else {
                        vdst = c_decoding_graph.get_vertex(BLUE_BOUNDARY_INDEX);
                    }
                } else {
                    vdst = c_decoding_graph.get_vertex(dst_id);
                }
                // We also need to handle the case when vthru is not in the restricted lattice.
                // However, we know that as color(src) = color(dst) != color(thru), there is a
                // vertex adjacent to all three of these vertices that is in the restricted lattice.
#ifdef DEBUG
                std::cout << "\tflag " << f << " is active: " << print_v(vsrc) 
                            << ", " << print_v(vthru) << ", " << print_v(vdst) << "\n";
#endif
                if (!c_decoding_graph[rc].contains(vthru)) {
                    bool found = false;
                    uint round_of_flag = src_id / detectors_per_round;
                    for (auto x : c_decoding_graph.get_neighbors(vsrc)) {
                        if (!is_colored_boundary(x)
                            && c_decoding_graph[rc].contains(x)
                            && c_decoding_graph.contains(x, vthru)
                            && c_decoding_graph.contains(x, vdst)
                            && (x->id / detectors_per_round) == round_of_flag)
                        {
                            vthru = x;
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
#ifdef DEBUG
                        std::cout << "\t\tfailed to find replacement vthru...\n";
#endif
                        continue;
                    }
#ifdef DEBUG
                    std::cout << "\t\treplaced vthru with " << print_v(vthru) << "\n";
#endif
                }
                DecodingGraph::flag_edge_t fe = std::make_tuple(vsrc, vthru, vdst);
                flag_edge_list.push_back(fe);
            }
            c_decoding_graph[rc].setup_flagged_decoding_graph(detector_list, flag_edge_list);
        }
        
        std::map<sptr<colored_vertex_t>, sptr<colored_vertex_t>> node_to_pref_boundary;

        for (uint i = 0; i < n; i++) {
            uint di = stab_events[i];
            auto vi = c_decoding_graph.get_vertex(di);
            for (uint j = i+1; j < n; j++) {
                uint dj = stab_events[j];
                sptr<colored_vertex_t> vj;
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

                    auto b1_error_data = get_error_chain_data(vi, vb1, rc);
                    auto b2_error_data = get_error_chain_data(vi, vb2, rc);
                    if (b1_error_data.weight < b2_error_data.weight) {
                        vj = vb1;
                    } else {
                        vj = vb2;
                    }
                    node_to_pref_boundary[vi] = vj;
                } else {
                    vj = c_decoding_graph.get_vertex(dj);
                }
                auto error_data = get_error_chain_data(vi, vj, rc);
                
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
            sptr<colored_vertex_t> bi;
            sptr<colored_vertex_t> bj;
            if (c_decoding_graph.are_matched_through_boundary(vi, vj, rcolors[k], &bi, &bj, flags_are_active)) {
                match_t m1 = make_match(vi->id, bi->id, rcolors[k]);
                match_t m2 = make_match(vj->id, bj->id, rcolors[k]);
                match_list.push_back(m1);
                match_list.push_back(m2);
            } else {
                match_t m = make_match(vi->id, vj->id, rcolors[k]);
                match_list.push_back(m);
            }
#ifdef DEBUG
            auto error_data = get_error_chain_data(vi, vj, rc);
            std::cout << "\t" << print_v(vi) << " <--> " << print_v(vj) << " (w = "
                << error_data.weight << ", thru b = " 
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

        sptr<e_t> e = std::make_shared<e_t>();
        e->src = (sptr<void>)v1;
        e->dst = (sptr<void>)v2;
        e->color = r;
        connectivity_graph.add_edge(e);
    }
    std::vector<cc_t> cc_list;
    // Compute connected components.
    auto rb = c_decoding_graph.get_vertex(RED_BOUNDARY_INDEX);  // Look from
    std::set<sptr<v_t>> skip_set;     // If we have a loop CC, then we need to skip certain searches.
    for (auto n : connectivity_graph.get_neighbors(rb)) {
        if (skip_set.count(n))  continue;
        sptr<e_t> rb_n_e = connectivity_graph.get_edge(rb, n);
        std::vector<match_t> match_list{make_match(RED_BOUNDARY_INDEX, n->id, rb_n_e->color)};

        std::map<sptr<v_t>, sptr<v_t>> prev;
        prev[n] = rb;
        auto curr = n;
        while (!is_colored_boundary(curr)) {
            // Get next neighbor.
            sptr<v_t> next = nullptr;
            for (auto w : connectivity_graph.get_neighbors(curr)) {
                if (w != prev[curr])    next = w;
            }
            if (next == nullptr)    break;
            prev[next] = curr;
            sptr<e_t> e = connectivity_graph.get_edge(curr, next);
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
    return cc_list;
}

void
RestrictionDecoder::insert_error_chain_into(
        std::set<sptr<colored_edge_t>>& edge_set,
        std::map<sptr<colored_edge_t>, uint>& incidence_map,
        std::string component_color,
        sptr<colored_vertex_t> src,
        sptr<colored_vertex_t> dst,
        std::string r_color)
{
#ifdef DEBUG
    std::cout << "\t" << print_v(src) << " <--> " << print_v(dst) << ":";
#endif
    auto error_data = get_error_chain_data(src, dst, r_color);
    auto error_chain = error_data.error_chain;
    for (uint j = 1; j < error_chain.size(); j++) {
        sptr<colored_vertex_t> u1 = std::static_pointer_cast<colored_vertex_t>(error_chain[j-1]);
        sptr<colored_vertex_t> u2 = std::static_pointer_cast<colored_vertex_t>(error_chain[j]);
        if (is_colored_boundary(u1) && is_colored_boundary(u2)) continue;
        if (u1->color != component_color && u2->color != component_color) continue;
        // Flatten the vertices.
        sptr<colored_vertex_t> fu1 = flatten(u1);
        sptr<colored_vertex_t> fu2 = flatten(u2);
        if (fu1 == fu2) continue;
        auto e = c_decoding_graph.get_edge(fu1, fu2);
        if (e == nullptr) {
            // Get path between the two vertices. This is likely some weird CX edge.
            insert_error_chain_into(edge_set, incidence_map, component_color, fu1, fu2, r_color);
        } else if (fu1->color == fu2->color) {
            // This is some weird CX edge. Get a common neighbor of the two vertices.
            sptr<colored_vertex_t> fu3;
            for (sptr<vertex_t> x : c_decoding_graph[r_color].get_common_neighbors(u1, u2)) {
                sptr<colored_vertex_t> _x = std::static_pointer_cast<colored_vertex_t>(x);
                sptr<colored_vertex_t> fx = flatten(_x);
                if (fx->color == fu1->color) continue;
                if (!c_decoding_graph.contains(fu1, fx) || !c_decoding_graph.contains(fu2, fx)) continue;
                fu3 = fx;
                break;
            }
            auto ex1 = c_decoding_graph.get_edge(fu1, fu3);
            auto ex2 = c_decoding_graph.get_edge(fu2, fu3);
#ifdef DEBUG
            std::cout << " " << print_v(fu1) << ", " << print_v(fu3) << ", " << print_v(fu2);
#endif
            edge_set.insert(ex1);
            edge_set.insert(ex2);
            incidence_map[ex1]++;
            incidence_map[ex2]++;
        } else {
#ifdef DEBUG
            std::cout << " " << print_v(fu1) << ", " << print_v(fu2);
#endif
            edge_set.insert(e);
            incidence_map[e]++;
        }
    }
#ifdef DEBUG
    std::cout << "\n";
#endif
}

sptr<colored_vertex_t>
RestrictionDecoder::flatten(sptr<colored_vertex_t> v) {
    if (is_colored_boundary(v)) return v;
    else {
        // In a circuit-level error model, flattening to the lowest level of detectors
        // (the first round) is not good as this ignores CNOT errors that occur in a prior
        // round. Thus, we instead flatten to the second round of detectors.
        uint64_t id = (v->id % detectors_per_round) + detectors_per_round;
        sptr<colored_vertex_t> fv = c_decoding_graph.get_vertex(id);
        return fv;
    }
}

std::set<face_t>
RestrictionDecoder::get_incident_faces(sptr<colored_vertex_t> v) {
    // We will be accessing this often. Just memoize.
    static std::map<sptr<colored_vertex_t>, std::set<face_t>> memo;
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
            if (!c_decoding_graph.contains_face(fc)) continue;
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

void
RestrictionDecoder::switch_out_boundaries(
        std::set<sptr<colored_edge_t>>& edge_set,
        std::map<sptr<colored_edge_t>, uint>& incidence_map) 
{
    std::set<sptr<colored_edge_t>> new_edges;
    for (auto it = edge_set.begin(); it != edge_set.end(); ) {
        sptr<colored_edge_t> e = *it;
        auto src = std::reinterpret_pointer_cast<colored_vertex_t>(e->src);
        auto dst = std::reinterpret_pointer_cast<colored_vertex_t>(e->dst);

        if (!is_colored_boundary(src) && !is_colored_boundary(dst)) {
            it++;
            continue;
        }
        
        sptr<colored_vertex_t> boundary, other;
        if (is_colored_boundary(src)) {
            boundary = src;
            other = dst;
        } else {
            boundary = dst;
            other = src;
        }

        uint64_t boundary_index;
        if (other->color == "r") {
            if (boundary->color == "g")  boundary_index = BLUE_BOUNDARY_INDEX;
            else                    boundary_index == GREEN_BOUNDARY_INDEX;
        } else if (other->color == "g") {
            if (boundary->color == "r")  boundary_index = BLUE_BOUNDARY_INDEX;
            else                    boundary_index == RED_BOUNDARY_INDEX;
        } else {
            if (boundary->color == "r")  boundary_index = GREEN_BOUNDARY_INDEX;
            else                    boundary_index == RED_BOUNDARY_INDEX;
        }
        auto new_boundary = c_decoding_graph.get_vertex(boundary_index);
        if (c_decoding_graph.contains(new_boundary, other)) {
            auto new_e = c_decoding_graph.get_edge(new_boundary, other);
            new_edges.insert(new_e);
            
            incidence_map[new_e] = incidence_map[e];
            incidence_map.erase(e);

            it = edge_set.erase(it);
        } else {
            it++;
        }
    }
    for (sptr<colored_edge_t> e : new_edges) edge_set.insert(e);
}

void
RestrictionDecoder::remap_boundary_edges(
        std::set<sptr<colored_edge_t>>& edge_set,
        std::map<sptr<colored_edge_t>, uint>& incidence_map) 
{
    std::set<sptr<colored_edge_t>> new_edges;
    for (auto it = edge_set.begin(); it != edge_set.end(); ) {
        sptr<colored_edge_t> e = *it;
        auto src = std::reinterpret_pointer_cast<colored_vertex_t>(e->src);
        auto dst = std::reinterpret_pointer_cast<colored_vertex_t>(e->dst);
        if (!is_colored_boundary(src) && !is_colored_boundary(dst)) {
            it++;
            continue;
        }
        sptr<colored_vertex_t> boundary, other;
        if (is_colored_boundary(src)) {
            boundary = src;
            other = dst;
        } else {
            boundary = dst;
            other = src;
        }
        // We want a common neighbor that is not the boundary, and has the same
        // color as the existing boundary.
        sptr<colored_vertex_t> common_non_boundary = nullptr;
        for (sptr<colored_vertex_t> x : c_decoding_graph.get_common_neighbors(boundary, other)) {
            if (is_colored_boundary(x)) continue;
            if (x->color != boundary->color) continue;
            sptr<colored_vertex_t> fx = flatten(x);
            if (!c_decoding_graph.contains(fx, other)) continue;
            common_non_boundary = fx;
            break;
        }

        if (common_non_boundary == nullptr) {
            it++;
            continue;
        }
        // Replace edge with new edge.
        sptr<colored_edge_t> new_e = c_decoding_graph.get_edge(other, common_non_boundary);
        if (new_e == nullptr) {
            it++;
        } else {
            new_edges.insert(new_e);

            incidence_map[new_e] = incidence_map[e];
            incidence_map.erase(e);

            it = edge_set.erase(it);
        }
    }
    for (auto e : new_edges) edge_set.insert(e);
}

}   // qontra
