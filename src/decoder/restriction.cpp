/*
 *  author: Suhas Vittal
 *  date:   19 July 2023
 * */

#include "decoder/restriction.h"

#include <assert.h>

#define DEBUG

namespace qontra {
namespace decoder {

using namespace graph;

RestrictionDecoder::RestrictionDecoder(const stim::Circuit& circuit)
    :Decoder(circuit, false),
    rlatt_dec(),
    to_rlatt(),
    from_rlatt(),
    color_map()
{
    // Restricted lattice order: RG (so blue is restricted), RB, GB
    const __COLOR colors[3] = { __COLOR::blue, __COLOR::green, __COLOR::red };
    // Create new Stim circuits for each restricted color.
    for (uint i = 0; i < 3; i++) {
        uint det_ctr = 0;
        uint subdet_ctr = 0;

        __COLOR restricted_color = colors[i];
        stim::Circuit sub_circuit;
        // Only keep detection events whose color is NOT blue.
        circuit.for_each_operation([&] (stim::Operation op) {
            std::string gate_name(op.gate->name);
            if (gate_name == "DETECTOR") {
                int color_id = (int) op.target_data.args[0];
                if (i2c(color_id) != restricted_color) {
                    sub_circuit.append_operation(op);

                    to_rlatt[det_ctr][i] = subdet_ctr;
                    from_rlatt[std::make_pair(subdet_ctr, i)] = det_ctr;
                    subdet_ctr++;
                }
                color_map[det_ctr++] = i2c(color_id);
            } else {
                sub_circuit.append_operation(op);
            }
        });
        rlatt_dec[i] = new MWPMDecoder(sub_circuit);
    }
}

Decoder::result_t
RestrictionDecoder::decode_error(const syndrome_t& syndrome) {
    const uint n_detectors = circuit.count_detectors();
    const uint n_observables = circuit.count_observables();

    timer.clk_start();
    std::vector<uint> detectors = get_nonzero_detectors(syndrome);

#ifdef DEBUG
    std::cout << "======================================\n";
    std::cout << "Detectors:";
    for (auto x : detectors) std::cout << " " << x << "(" << c2i(color_map[x]) << ")";
    std::cout << "\n";
#endif

    // Get matching results for each restricted lattice.
    auto r_lattRG = decode_restricted_lattice(detectors, __COLOR::blue);
    auto r_lattRB = decode_restricted_lattice(detectors, __COLOR::green);
    auto r_lattGB = decode_restricted_lattice(detectors, __COLOR::red);

    typedef std::pair<uint, __COLOR>    cdet_t;
    // Compute the connections for each detector. Also, we need to track
    // their color for cases such as the boundary (where the color is
    // pre-determined).
    // 
    // Create a graph to track the connections and then later compute the
    // connected components.
    struct cv_t : base::vertex_t {
        cdet_t  detector;
    };
    struct ce_t : base::edge_t {
        __COLOR color;
    };
    typedef Graph<cv_t, ce_t>   ConnGraph;

    ConnGraph connection_graph;
    Decoder::result_t* res_array[] = { &r_lattRG, &r_lattRB, &r_lattGB };

    const __COLOR restricted_colors[] = {__COLOR::blue, __COLOR::green, __COLOR::red};
    for (uint i = 0; i < 3; i++) {
        auto res_p = res_array[i];
        for (auto aa : res_p->error_assignments) {
            uint d1 = std::get<0>(aa);
            uint d2 = std::get<1>(aa);

            if (d1 == BOUNDARY_INDEX)   std::swap(d1, d2);
            __COLOR c1 = color_map[d1];
            // Compute d2's color. If d2 is the boundary, then we need
            // to use d1's color to do so.
            __COLOR c2;
            if (d2 == BOUNDARY_INDEX) {
                uint rlatt_d1 = to_rlatt[d1][i];
                // Get decoding graph vertices.
                DecodingGraph& gr = rlatt_dec[i]->decoding_graph;
                auto vd1 = gr.get_vertex(rlatt_d1);
                auto vd2 = gr.get_vertex(BOUNDARY_INDEX);
                auto ec_data = gr.get_error_chain_data(vd1, vd2);
                // If the error chain length is even, then d2 is the same
                // color as d1. Otherwise, it is the opposite color.
                uint len = ec_data.chain_length;
                if (len & 0x1) {
                    c2 = get_remaining_color(c1, restricted_colors[i]);
                } else {
                    c2 = c1;
                }
            } else {
                c2 = color_map[d2];
            }
            cdet_t cd1 = std::make_pair(d1, c1);
            cdet_t cd2 = std::make_pair(d2, c2);
            // Create an edge in the connection graph
            // Create any vertices if they do not exist.
            cv_t* v1;
            cv_t* v2;
            uint64_t id1 = cdet_to_id(cd1);
            uint64_t id2 = cdet_to_id(cd2);
            if ((v1=connection_graph.get_vertex(id1)) == nullptr) {
                v1 = new cv_t;
                v1->id = id1;
                v1->detector = cd1;
                connection_graph.add_vertex(v1);
            }
            if ((v2=connection_graph.get_vertex(id2)) == nullptr) {
                v2 = new cv_t;
                v2->id = id2;
                v2->detector = cd2;
                connection_graph.add_vertex(v2);
            }
            ce_t* e = new ce_t;
            e->src = v1;
            e->dst = v2;
            e->is_undirected = true;
            e->color = restricted_colors[i];
            connection_graph.add_edge(e);
        }
    }
    // Here, we will partition the syndrome into two parts:
    //  (1) Syndrome bits not in a connected component with the red boundary.
    //  (2) Everything else.
    cdet_t red_boundary = std::make_pair(BOUNDARY_INDEX, __COLOR::red);
    auto rbv = connection_graph.get_vertex(cdet_to_id(red_boundary));

    std::set<cdet_t> in_connected_components{red_boundary};

    std::vector<component_t> connected_components;
    if (rbv != nullptr) {
        // Perform a DFS to get paths from the red boundary to any other boundary.
        //
        // Needs to be a special DFS, so we aren't using the xfs function.
        std::deque<cv_t*> dfs;
        dfs.push_back(rbv);

        std::set<cv_t*> visited;
        std::map<cv_t*, cv_t*> prev;
        while (dfs.size()) {
            auto v = dfs.back();
            dfs.pop_back();
            if (visited.count(v))   continue;
            visited.insert(v);
            for (auto w : connection_graph.get_neighbors(v)) {
                if (w == prev[v])   continue;
                cdet_t wd = w->detector;
                if (wd.first == BOUNDARY_INDEX) {
                    // This is one connected component.
                    std::vector<match_t> cc;

                    auto e_vw = connection_graph.get_edge(v, w);
                    cc.push_back(std::make_tuple(wd, v->detector, e_vw->color));

                    auto curr = v;
                    while (curr != rbv) {
                        in_connected_components.insert(curr->detector);
                        auto p = prev[curr];
                        auto e = connection_graph.get_edge(curr, p);
                        cc.push_back(std::make_tuple(p->detector, curr->detector, e->color));
                        curr = p;
                    }
                    __COLOR cc_color = get_remaining_color(__COLOR::red, wd.second);
                    connected_components.push_back(std::make_pair(cc, cc_color));
                } else if (!visited.count(w)) {
                    prev[w] = v;
                    dfs.push_back(w);
                }
            }
        }
    }
    // Add everything not in a connected components to its own component.
    std::vector<match_t> not_cc;
    std::set<cdet_t> visited;
    for (uint d : detectors) {
        __COLOR c = color_map[d];
        cdet_t det = std::make_pair(d, c);
        if (in_connected_components.count(det))    continue;

        auto v = connection_graph.get_vertex(cdet_to_id(det));
        for (auto w : connection_graph.get_neighbors(v)) {
            if (visited.count(w->detector))   continue;        
            if (in_connected_components.count(w->detector)) continue;
            auto e = connection_graph.get_edge(v, w);
            not_cc.push_back(std::make_tuple(det, w->detector, e->color));
        }
        visited.insert(det);
    }
    connected_components.push_back(std::make_pair(not_cc, __COLOR::red));
    stim::simd_bits corr(n_observables);
    corr.clear();

#ifdef DEBUG
    std::cout << "Connected Components:\n";
#endif
    for (auto& pair : connected_components) {
        auto cc = pair.first;
        auto cc_color = pair.second;

        std::set<cdetpair_t> in_face;

#ifdef DEBUG
        std::cout << "\tColor " << c2i(cc_color) << ":\n";
#endif
        
        auto flipped_edges = get_all_edges_in_component(cc);
        std::set<cdet_t> vertices_in_flips; // of the same color as the component.
        std::map<cdet_t, cdet_t> starting_neighbor;
#ifdef DEBUG
        std::cout << "\t\tEdges:\n";
#endif
        for (cdetpair_t edge : flipped_edges) { 
            cdet_t cdx = edge.first;
            cdet_t cdy = edge.second;
#ifdef DEBUG
            std::cout << "\t\t" << cdx.first << "(" << c2i(cdx.second) << ")"
                        << ", " << cdy.first << "(" << c2i(cdy.second) << ")\n";
#endif
            if (cdx.second == cdy.second)   continue;
            if (cdx.second == cc_color) {
                vertices_in_flips.insert(cdx);
                starting_neighbor[cdx] = cdy;
            } 
            if (cdy.second == cc_color) {
                vertices_in_flips.insert(cdy);
                starting_neighbor[cdy] = cdx;
            }
        }
        // Now, we will apply corrections on faces around each vertex.
        for (cdet_t cdx : vertices_in_flips) {
#ifdef DEBUG
            std::cout << "\t\tChecking vertex " << cdx.first << "(" << c2i(cdx.second) << "):\n";
#endif
            auto neighbors = get_neighbors(cdx);
#ifdef DEBUG
            std::cout << "\t\t\tneighbors:";
            for (auto x : neighbors) std::cout << " " << x.first << "(" << c2i(x.second) << ")";
            std::cout << "\n";
#endif
            auto cdy = starting_neighbor[cdx];
            while (neighbors.size()) {
                // Get common neighbor of cdx and cdy.
                neighbors.erase(cdy);
                auto n = get_neighbors(cdy);
                std::set<cdet_t> common;
                std::set_intersection(neighbors.begin(), neighbors.end(),
                                        n.begin(), n.end(),
                                        std::inserter(common, common.begin()));
                if (common.empty()) break;
                auto cdz = *common.begin();
                auto e_xy = std::make_pair(cdx, cdy);
                auto e_xz = std::make_pair(cdx, cdz);
                auto e_yz = std::make_pair(cdy, cdz);

                auto e_yx = std::make_pair(cdy, cdx);
                auto e_zx = std::make_pair(cdz, cdx);
                auto e_zy = std::make_pair(cdz, cdy);

                if (in_face.count(e_xy) || in_face.count(e_yx)
                    || in_face.count(e_xz) || in_face.count(e_zx)
                    || in_face.count(e_yz) || in_face.count(e_zy))
                {
                    break;
                }
                corr ^= get_correction_for_face(cdx, cdy, cdz);
                if (flipped_edges.count(e_xy) || flipped_edges.count(e_yx)) {
                    in_face.insert(e_xy);
                    in_face.insert(e_yx);
                }
                if (flipped_edges.count(e_xz) || flipped_edges.count(e_zx)) {
                    in_face.insert(e_xz);
                    in_face.insert(e_zx);
                }
                if (flipped_edges.count(e_yz) || flipped_edges.count(e_zy)) {
                    in_face.insert(e_yz);
                    in_face.insert(e_zy);
                }
                cdy = cdz;
            }
        }
    }

#ifdef DEBUG
    std::cout << "\tis error : " << is_error(corr, syndrome) << "\n";
#endif
    fp_t t = (fp_t)timer.clk_end();
    return (Decoder::result_t) { t, corr, is_error(corr, syndrome) };
}

Decoder::result_t
RestrictionDecoder::decode_restricted_lattice(
                            const std::vector<uint>& detectors,
                            __COLOR restricted_color)
{
    const uint i = 2 - c2i(restricted_color);
    MWPMDecoder* dec = rlatt_dec[i];
    const uint n_detectors = dec->get_circuit().count_detectors();
    const uint n_observables = dec->get_circuit().count_observables();

    syndrome_t restricted_syndrome(n_detectors+n_observables);
    restricted_syndrome.clear();
    for (uint d : detectors) {
        if (color_map[d] == restricted_color)   continue;
        uint dd = to_rlatt[d][i];
        restricted_syndrome[dd] = 1;
    }
    auto res = dec->decode_error(restricted_syndrome);
    // Convert the error assignments to the original circuit's detectors.
#ifdef DEBUG
    std::cout << "\tAssignments on lattice restricting color " 
        << c2i(restricted_color) << ":\n";
#endif
    // First, filter out any assignments that go through the boundary.
    std::vector<Decoder::assign_t>  new_assignments;
    DecodingGraph& gr = dec->decoding_graph;
    for (auto aa : res.error_assignments) {
        uint d1 = std::get<0>(aa);
        uint d2 = std::get<1>(aa);
        if (d1 != BOUNDARY_INDEX && d2 != BOUNDARY_INDEX) {
            auto v1 = gr.get_vertex(d1);
            auto v2 = gr.get_vertex(d2);
            if (gr.get_error_chain_data(v1, v2).error_chain_runs_through_boundary) {
                new_assignments.push_back(std::make_tuple(d1, BOUNDARY_INDEX, std::get<2>(aa)));
                new_assignments.push_back(std::make_tuple(d2, BOUNDARY_INDEX, std::get<2>(aa)));
            } else {
                new_assignments.push_back(aa);
            }
        } else {
            new_assignments.push_back(aa);
        }
    }
    for (auto& aa : new_assignments) {
        uint& d1 = std::get<0>(aa);
        uint& d2 = std::get<1>(aa);
        if (d1 != BOUNDARY_INDEX)   d1 = from_rlatt[std::make_pair(d1, i)];
        if (d2 != BOUNDARY_INDEX)   d2 = from_rlatt[std::make_pair(d2, i)];
#ifdef DEBUG
        std::cout << "\t\t" << d1 << " <----> " << d2 << "\n";
#endif
    }
    res.error_assignments = new_assignments;
    return res;
}

std::set<cdetpair_t>
RestrictionDecoder::get_all_edges_in_component(const std::vector<match_t>& comp) {
    std::set<cdetpair_t> edges;
    for (match_t mate : comp) {
        cdet_t cdx = std::get<0>(mate);
        cdet_t cdy = std::get<1>(mate);
        __COLOR restricted_color = std::get<2>(mate);

        uint dec_index = 2 - c2i(restricted_color);
        DecodingGraph& gr = rlatt_dec[dec_index]->decoding_graph;
        uint64_t ddx = cdx.first == BOUNDARY_INDEX ? cdx.first : to_rlatt[cdx.first][dec_index];
        uint64_t ddy = cdy.first == BOUNDARY_INDEX ? cdy.first : to_rlatt[cdy.first][dec_index];

        auto vx = gr.get_vertex(ddx);
        auto vy = gr.get_vertex(ddy);
        auto error_chain = gr.get_error_chain_data(vx, vy).error_chain;
        for (uint i = 1; i < error_chain.size(); i++) {
            auto v1 = error_chain[i-1];
            auto v2 = error_chain[i];

            uint64_t d1 = v1->id;
            uint64_t d2 = v2->id;
            if (d1 == BOUNDARY_INDEX)   std::swap(d1, d2);
            d1 = from_rlatt[std::make_pair(d1, dec_index)];
            __COLOR c1 = color_map[d1];
            cdet_t cd1 = std::make_pair(d1, c1);
            __COLOR c2;
            if (d2 == BOUNDARY_INDEX) {
                c2 = get_remaining_color(c1, restricted_color);
            } else {
                d2 = from_rlatt[std::make_pair(d2, dec_index)];
                c2 = color_map[d2];
            }
            cdet_t cd2 = std::make_pair(d2, c2);
            edges.insert(std::make_pair(cd1, cd2));
        }
    }
    return edges;
}

std::set<cdet_t>
RestrictionDecoder::get_common_neighbors(cdet_t cdx, cdet_t cdy) {
    auto adjx = get_neighbors(cdx);
    auto adjy = get_neighbors(cdy);

    std::set<cdet_t> common;
    std::set_intersection(adjx.begin(), adjx.end(), adjy.begin(), adjy.end(),
                            std::inserter(common, common.begin()));
    // Filter out all detectors with the same color as either cdx or cdy
    for (auto it = common.begin(); it != common.end(); ) {
        if (it->second == cdx.second || it->second == cdy.second)   it = common.erase(it);
        else                                                        it++;
    }
    return common;
}

std::set<cdet_t>
RestrictionDecoder::get_neighbors(cdet_t cd) {
    std::set<cdet_t> neighbors;
    for (uint i = 0; i < 3; i++) {
        if (c2i(cd.second) == i)    continue;
        __COLOR restricted_color = i2c(i);
        uint dec_index = 2 - i;
        DecodingGraph& gr = rlatt_dec[dec_index]->decoding_graph;

        uint64_t dv = cd.first == BOUNDARY_INDEX ? BOUNDARY_INDEX : to_rlatt[cd.first][dec_index];
        auto v = gr.get_vertex(dv);
        for (auto w : gr.get_neighbors(v)) {
            uint64_t dw = w->id;
            __COLOR cw;
            if (dw == BOUNDARY_INDEX) {
                cw = get_remaining_color(cd.second, restricted_color);
            } else {
                dw = from_rlatt[std::make_pair(dw, dec_index)];
                cw = color_map[dw];
            }
            if (cd.second == cw)    continue;
            neighbors.insert(std::make_pair(dw, cw));
        }
        if (cd.first == BOUNDARY_INDEX) {
            __COLOR cb = get_remaining_color(cd.second, restricted_color);
            neighbors.insert(std::make_pair(BOUNDARY_INDEX, cb));
        }
    }
    return neighbors;
}

stim::simd_bits
RestrictionDecoder::get_correction_for_face(cdet_t x, cdet_t y, cdet_t z) {
    const std::vector<cdet_t> detectors{x, y, z};

    stim::simd_bits corr(circuit.count_observables());
    corr.clear();
#ifdef DEBUG
    std::cout << "\t\t\tApplying corrections on face " 
                << x.first << "(" << c2i(x.second) << "), "
                << y.first << "(" << c2i(y.second) << "), "
                << z.first << "(" << c2i(z.second) << "):";
#endif
    for (uint i = 0; i < detectors.size(); i++) {
        cdet_t cd1 = detectors[i];
        for (uint j = i+1; j < detectors.size(); j++) {
            cdet_t cd2 = detectors[j];
            __COLOR restricted_color = get_remaining_color(cd1.second, cd2.second);
            uint dec_index = 2 - c2i(restricted_color);
            DecodingGraph& gr = rlatt_dec[dec_index]->decoding_graph;

            uint64_t dd1 = cd1.first == BOUNDARY_INDEX 
                                ? BOUNDARY_INDEX : to_rlatt[cd1.first][dec_index];
            uint64_t dd2 = cd2.first == BOUNDARY_INDEX 
                                ? BOUNDARY_INDEX : to_rlatt[cd2.first][dec_index];
            auto v1 = gr.get_vertex(dd1);
            auto v2 = gr.get_vertex(dd2);
            for (auto fr : gr.get_error_chain_data(v1, v2).frame_changes) {
                corr[fr] ^= 1;
#ifdef DEBUG
                std::cout << " FR" << fr;
#endif
            }
        }
    }
#ifdef DEBUG
    std::cout << "\n";
#endif
    return corr;
}

}   // decoder
}   // qontra
