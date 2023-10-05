/*
 *  author: Suhas Vittal
 *  date:   19 July 2023
 * */

#include "decoder/restriction.h"

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
        std::vector<match_t> matches = std::get<3>(comp);
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
        in_cc_array.push_back(std::make_pair(es, std::get<4>(comp)));
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
            if (!all_in_cc.count(e)) es.insert(e);
        }
    }
    
    std::set<colored_vertex_t*> out_of_cc_incident = get_incident_vertices(out_of_cc, "r");
    std::set<colored_vertex_t*> in_cc_incident;
    for (edge_set_t& es : in_cc_array) {
        auto incident = get_incident_vertices(es.first, es.second);
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
        assert(nf <= 16);   // Exit, because otherwise we will take too long.
        uint64_t enf = 1L << nf;

        uint64_t best_intersect_with_cc = 0;
        uint64_t best_intersect_with_no_cc = 0;
        std::set<colored_edge_t> best_boundary;
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

                    local_corr ^= c_decoding_graph.get_correction_for_face(f);
                }
                j >>= 1;
                it++;
            }
            // Check how much the edges on the boundary of the
            // faces intersect with in_cc or out_of_cc.
            std::set<colored_edge_t*> int_with_cc, int_with_no_cc;
            std::set_intersection(
                    f_boundary.begin(), f_boundary.end(),
                    in_cc.begin(), in_cc.end(),
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
            in_cc.erase(e);
            out_of_cc.erase(e);
        }
    }

    return (Decoder::result_t) {
        0.0,
        corr,
        is_error(syndrome, corr)
    };
}

std::vector<match_t>
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
            occurrence_by_color[0]++;
        } else if (v->color == "g") {
            restricted_syndromes[0].push_back(d);
            restricted_syndromes[2].push_back(d);
            occurrence_by_color[1]++;
        } else {
            restricted_syndromes[1].push_back(d);
            restricted_syndromes[2].push_back(d);
            occurrence_by_color[2]++;
        }
    }

    // Add boundaries if necessary.
    for (std::vector<uint>& rs : restricted_syndromes) {
        if (rs.size() & 1)  rs.push_back(BOUNDARY_INDEX);
    }

    // Perform MWPM. As this is a bit specialized, we just implement by hand.
    const std::string rcolors[] = { "rg", "rb", "gb" };
    for (uint k = 0; k < 3; k++) {
        const uint n = restricted_syndromes[k].size();
        const uint m = (n*(n+1))/2;
        PerfectMatching pm(n, m);
        for (uint i = 0; i < n; i++) {
            uint di = restricted_syndromes[k][i];
            auto vi = c_decoding_graph.get_vertex(di);
            for (uint j = i+1; j < n; j++) {
                uint dj = restricted_syndromes[k][j];
                auto vj = c_decoding_graph.get_vertex(dj);
                auto error_data = c_decoding_graph[rcolors[k]].get_error_chain_data(vi, vj);
                
                wgt_t edge_weight;
                if (error_data.weight > 1000.0) edge_weight = 1000000000;
                else                            edge_weight = (wgt_t) 1000.0*error_data.weight;
                pm.AddEdge(i, j, edge_weight);
            }
        }
        pm.Solve();
        for (uint i = 0; i < n; i++) {
            uint j = pm.GetMatch(i);
            uint di = restricted_syndromes[k][i],
                dj = restricted_syndromes[k][j];
        }
    }
}



}   // qontra
