/*
 *  author: Suhas Vittal
 *  date:   19 July 2023
 * */

#include "decoder/restriction.h"

#define DEBUG

namespace qontra {

using namespace graph;
using namespace decoding;

template <class T> void
xor_entry_into(T x, std::set<T>& s) {
    if (s.count(x)) s.erase(x);
    else            s.insert(x);
}

Decoder::result_t
RestrictionDecoder::decode_error(const syndrome_t& syndrome) {
#ifdef DEBUG
    std::cout << "===================================\n";
#endif

    std::vector<uint> detectors = get_nonzero_detectors(syndrome);
    std::vector<match_t> matchings = blossom_subroutine(detectors);
    std::vector<cc_t> comps = compute_connected_components(matchings);
    // Now, divide the edges between those in the connected components
    // and those not in the components.
    typedef std::pair<std::set<colored_edge_t*>, std::string> edge_set_t;

    std::set<colored_edge_t*> out_of_cc;
    std::vector<edge_set_t> in_cc_array;
    std::set<colored_edge_t*> all_in_cc;

    // Compute in cc.
    for (cc_t& comp : comps) {
        std::vector<match_t> matches = std::get<2>(comp);
        std::set<colored_edge_t*> es;
        for (match_t& m : matches) {
            auto v = c_decoding_graph.get_vertex(std::get<0>(m));
            auto w = c_decoding_graph.get_vertex(std::get<1>(m));
            std::string r = std::get<2>(m);
            auto error_data = c_decoding_graph[r].get_error_chain_data(v, w);
            auto error_chain = error_data.error_chain;
            for (uint j = 1; j < error_chain.size(); j++) {
                colored_vertex_t* u1 = (colored_vertex_t*)error_chain[j-1];
                colored_vertex_t* u2 = (colored_vertex_t*)error_chain[j];
                auto e = c_decoding_graph.get_edge(u1, u2);
                es.insert(e);
                all_in_cc.insert(e);
            }
        }
        in_cc_array.push_back(std::make_pair(es, std::get<3>(comp)));
    }
    // Now compute out of cc, which contains anything not in cc.
    for (match_t& m : matchings) {
        colored_vertex_t* v = c_decoding_graph.get_vertex(std::get<0>(m));
        colored_vertex_t* w = c_decoding_graph.get_vertex(std::get<1>(m));
        std::string r = std::get<2>(m);
        auto error_data = c_decoding_graph[r].get_error_chain_data(v, w);
        auto error_chain = error_data.error_chain;
        for (uint j = 1; j < error_chain.size(); j++) {
            colored_vertex_t* u1 = (colored_vertex_t*)error_chain[j-1];
            colored_vertex_t* u2 = (colored_vertex_t*)error_chain[j];
            auto e = c_decoding_graph.get_edge(u1, u2);
            if (!all_in_cc.count(e)) out_of_cc.insert(e);
        }
    }
    
    std::set<colored_vertex_t*> out_of_cc_incident = c_decoding_graph.get_all_incident_vertices(out_of_cc, "r");
    std::set<colored_vertex_t*> in_cc_incident;
    for (edge_set_t& es : in_cc_array) {
        auto incident = c_decoding_graph.get_all_incident_vertices(es.first, es.second);
        for (auto x : incident) in_cc_incident.insert(x);
    }
    std::set<colored_vertex_t*> all_incident;
    for (auto v : out_of_cc_incident)   all_incident.insert(v);
    for (auto v : in_cc_incident)       all_incident.insert(v);

    // Finally, compute the correction.
    const uint obs = circuit.count_observables();
    stim::simd_bits corr(obs);
    corr.clear();

    for (auto v : all_incident) {
        // Can only match one of out_of_cc or in_cc.
        std::set<face_t> incident_faces = c_decoding_graph.get_all_incident_faces(v);
        uint64_t nf = incident_faces.size();
        if (nf > 16) {
            std::cout << "nf = " << nf << "\n";
            assert(nf <= 16);   // Exit, because otherwise we will take too long.
        }
        uint64_t enf = 1L << nf;

        uint64_t best_intersect_with_cc = 0;
        uint64_t best_intersect_with_no_cc = 0;
        std::set<colored_edge_t*> best_boundary;
        stim::simd_bits best_corr(obs);
        for (uint64_t i = 0; i < enf; i++) {
            // Interpret the bits of i as the faces we will examine.
            std::set<colored_edge_t*> f_boundary;
            uint64_t j = i;
            auto it = incident_faces.begin();

            stim::simd_bits local_corr(obs);
            local_corr.clear();
            while (j) {
                if (j & 1) {
                    face_t f = *it;
                    colored_vertex_t* v1 = std::get<0>(f);
                    colored_vertex_t* v2 = std::get<1>(f);
                    colored_vertex_t* v3 = std::get<2>(f);

                    auto e12 = c_decoding_graph.get_edge(v1, v2);
                    auto e13 = c_decoding_graph.get_edge(v1, v3);
                    auto e23 = c_decoding_graph.get_edge(v2, v3);
                    xor_entry_into(e12, f_boundary);
                    xor_entry_into(e13, f_boundary);
                    xor_entry_into(e23, f_boundary);

                    local_corr ^= get_correction_for_face(f);
                }
                j >>= 1;
                it++;
            }
            // Check how much the edges on the boundary of the
            // faces intersect with in_cc or out_of_cc.
            std::set<colored_edge_t*> int_with_cc, int_with_no_cc;
            std::set_intersection(
                    f_boundary.begin(), f_boundary.end(),
                    all_in_cc.begin(), all_in_cc.end(),
                    std::inserter(int_with_cc, int_with_cc.begin()));
            std::set_intersection(
                    f_boundary.begin(), f_boundary.end(),
                    out_of_cc.begin(), out_of_cc.end(),
                    std::inserter(int_with_no_cc, int_with_no_cc.begin()));
            if ((best_intersect_with_cc >= best_intersect_with_no_cc
                    && int_with_cc.size() > best_intersect_with_cc)
                || (best_intersect_with_no_cc > best_intersect_with_cc
                    && int_with_no_cc.size() > best_intersect_with_no_cc))
            {
                best_intersect_with_cc = int_with_cc.size();
                best_intersect_with_no_cc = int_with_no_cc.size();
                best_boundary = f_boundary;
                best_corr = local_corr;
            }
        }
        // Commit the correction for the boundary and erase the edges.
        corr ^= best_corr;
        for (auto e : best_boundary) {
            all_in_cc.erase(e);
            out_of_cc.erase(e);
        }
    }

#ifdef DEBUG
    std::cout << "is error : " << is_error(syndrome, corr) << "\n";
#endif

    return (Decoder::result_t) {
        0.0,
        corr,
        is_error(syndrome, corr)
    };
}

std::vector<RestrictionDecoder::match_t>
RestrictionDecoder::blossom_subroutine(const std::vector<uint>& detectors) {
    // Partition the detectors into syndromes for each restricted lattice and
    // compute the MWPM.
    std::array<std::vector<uint>, 3> restricted_syndromes;
    restricted_syndromes.fill(std::vector<uint>());
    for (uint d : detectors) {
        auto v = c_decoding_graph.get_vertex(d);
        if (v->color == "r") {
            restricted_syndromes[0].push_back(d);
            restricted_syndromes[1].push_back(d);
        } else if (v->color == "g") {
            restricted_syndromes[0].push_back(d);
            restricted_syndromes[2].push_back(d);
        } else {
            restricted_syndromes[1].push_back(d);
            restricted_syndromes[2].push_back(d);
        }
    }

    // Add boundaries if necessary.
    for (std::vector<uint>& rs : restricted_syndromes) {
        if (rs.size() & 1)  rs.push_back(BOUNDARY_INDEX);
    }

    // Perform MWPM. As this is a bit specialized, we just implement by hand.
    const std::string rcolors[] = { "rg", "rb", "gb" };
    std::vector<match_t> match_list;
    for (uint k = 0; k < 3; k++) {
        const uint n = restricted_syndromes[k].size();
        const uint m = (n*(n+1))/2;
        PerfectMatching pm(n, m);
        pm.options.verbose = false;
        
        std::map<colored_vertex_t*, colored_vertex_t*> node_to_pref_boundary;

        for (uint i = 0; i < n; i++) {
            uint di = restricted_syndromes[k][i];
            auto vi = c_decoding_graph.get_vertex(di);
            for (uint j = i+1; j < n; j++) {
                uint dj = restricted_syndromes[k][j];
                colored_vertex_t* vj;
                if (dj == BOUNDARY_INDEX) {
                    // We need to get the correct boundary, which is the one
                    // closer to vi.
                    uint64_t b1, b2;
                    if (rcolors[k] == "r")  b1 = RED_BOUNDARY_INDEX;
                    else                    b2 = GREEN_BOUNDARY_INDEX;

                    if (rcolors[k] == "g")  b1 = GREEN_BOUNDARY_INDEX;
                    else                    b2 = BLUE_BOUNDARY_INDEX;

                    auto vb1 = c_decoding_graph.get_vertex(b1);
                    auto vb2 = c_decoding_graph.get_vertex(b2);

                    auto b1_error_data = c_decoding_graph[rcolors[k]].get_error_chain_data(vi, vb1);
                    auto b2_error_data = c_decoding_graph[rcolors[k]].get_error_chain_data(vi, vb2);
                    if (b1_error_data.weight < b2_error_data.weight) {
                        vj = vb1;
                    } else {
                        vj = vb2;
                    }
                    node_to_pref_boundary[vi] = vj;
                } else {
                    vj = c_decoding_graph.get_vertex(dj);
                }
                auto error_data = c_decoding_graph[rcolors[k]].get_error_chain_data(vi, vj);
                
                uint32_t edge_weight;
                if (error_data.weight > 1000.0) edge_weight = 1000000000;
                else                            edge_weight = (uint32_t) 1000.0*error_data.weight;
                pm.AddEdge(i, j, edge_weight);
            }
        }
        pm.Solve();

#ifdef DEBUG
        std::cout << "For restricted lattice " << k << ":\n";
#endif
        for (uint i = 0; i < n; i++) {
            uint j = pm.GetMatch(i);
            if (i >= j) continue;
            uint di = restricted_syndromes[k][i],
                dj = restricted_syndromes[k][j];
            auto vi = c_decoding_graph.get_vertex(di);
            auto vj = dj == BOUNDARY_INDEX ? node_to_pref_boundary[vi] : c_decoding_graph.get_vertex(dj);
            match_t m = std::make_tuple(vi->id, vj->id, rcolors[k]);
#ifdef DEBUG
            auto error_data = c_decoding_graph[rcolors[k]].get_error_chain_data(vi, vj);
            std::cout << "\t" << vi->id << "[" << vi->color << "] <--> " 
                << vj->id << "[" << vj->color << "] (thru b = " 
                << error_data.error_chain_runs_through_boundary << ")\n";
#endif
            match_list.push_back(m);
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

        colored_vertex_t* b1;
        colored_vertex_t* b2;
        if (c_decoding_graph.are_matched_through_boundary(v1, v2, r, &b1, &b2)) {
            e_t* e1 = new e_t;
            e_t* e2 = new e_t;
            e1->src = (void*)v1;
            e1->dst = (void*)b1;
            e1->color = r;
            e2->src = (void*)v2;
            e2->dst = (void*)b2;
            e2->color = r;
            connectivity_graph.add_edge(e1);
            connectivity_graph.add_edge(e2);
        } else {
            e_t* e = new e_t;
            e->src = (void*)v1;
            e->dst = (void*)v2;
            e->color = r;
            connectivity_graph.add_edge(e);
        }
    }
    std::vector<cc_t> cc_list;
    // Compute connected components.
    auto rb = c_decoding_graph.get_vertex(RED_BOUNDARY_INDEX);  // Look from
    for (auto n : connectivity_graph.get_neighbors(rb)) {
        e_t* rb_n_e = connectivity_graph.get_edge(rb, n);
        std::vector<match_t> match_list{std::make_tuple(RED_BOUNDARY_INDEX, n->id, rb_n_e->color)};

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
            match_list.push_back(std::make_tuple(curr->id, next->id, e->color));
            curr = next;
        }
        // Only proceeded if we ended up at a boundary.
        if (is_colored_boundary(curr)) {
            uint64_t b1 = RED_BOUNDARY_INDEX;
            uint64_t b2 = curr->id;

            std::string cc_color;
            if (b2 == RED_BOUNDARY_INDEX || b2 == BLUE_BOUNDARY_INDEX) {
                cc_color = "g";
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

stim::simd_bits
RestrictionDecoder::get_correction_for_face(face_t fc) {
    // Here, we simply AND all the correction bits together.
    auto v1 = std::get<0>(fc);
    auto v2 = std::get<1>(fc);
    auto v3 = std::get<2>(fc);
    
    const uint obs = circuit.count_observables();
    stim::simd_bits corr(obs);
    corr.clear();
    corr.invert_bits(); // Should be all ones now.

    colored_vertex_t* varray[] = {v1, v2, v3};
    for (uint i = 0; i < 3; i++) {
        auto vi = varray[i];
        for (uint j = i+1; j < 3; j++) {
            auto vj = varray[j];
            // An exception to this strategy is when vi and vj are both
            // boundaries. Here, they are merely used to complete the face,
            // and do not flip any observables. Skip them in this situation.
            if (is_colored_boundary(vi) && is_colored_boundary(vj)) continue;

            std::string r = vi->color + vj->color;
            auto error_data = c_decoding_graph[r].get_error_chain_data(vi, vj);
            stim::simd_bits local_corr(obs);
            local_corr.clear();
            for (uint fr : error_data.frame_changes) {
                local_corr[fr] = 1;
            }
            corr &= local_corr;
        }
    }
    return corr;
}

}   // qontra
