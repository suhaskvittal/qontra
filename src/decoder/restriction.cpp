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
                if (colors[color_id] != restricted_color) {
                    sub_circuit.append_operation(op);

                    to_rlatt[det_ctr][i] = subdet_ctr;
                    from_rlatt[std::make_pair(subdet_ctr, i)] = det_ctr;
                    subdet_ctr++;
                }
                color_map[det_ctr++] = colors[color_id];
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
    const __COLOR rlatt_c1_array[] = { __COLOR::red, __COLOR::red, __COLOR::green };
    const __COLOR rlatt_c2_array[] = { __COLOR::green, __COLOR::blue, __COLOR::blue };
    for (uint i = 0; i < 3; i++) {
        auto res_p = res_array[i];
        __COLOR first_possible_color = rlatt_c1_array[i];
        __COLOR second_possible_color = rlatt_c2_array[i];

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
                    if (c1 == first_possible_color) c2 = second_possible_color;
                    else                            c2 = first_possible_color;
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
    // Then, we will "jam" matchings in each partition separately. Jamming is essentially
    // when a detector has matchings in two different colors (i.e. G and B). Then, we ignore
    // one of the incident frame changes as this implies that we are at a face.
    typedef std::tuple<cdet_t, cdet_t, __COLOR>         match_t;
    typedef std::pair<std::vector<match_t>, __COLOR>    jamming_set_t;
    
    cdet_t red_boundary = std::make_pair(BOUNDARY_INDEX, __COLOR::red);
    auto rbv = connection_graph.get_vertex(cdet_to_id(red_boundary));

    std::set<cdet_t> in_connected_components{red_boundary};

    std::vector<jamming_set_t> jamming_sets;
    if (rbv != nullptr) {
        // Perform a DFS to get paths from the red boundary to any other boundary.
        //
        // Needs to be a special DFS, so we aren't using the xfs function.
        std::deque<cv_t*> dfs;
        dfs.push_back(rbv);

        std::set<cv_t*> visited;
        std::map<cv_t*, cv_t*> prev;
#ifdef DEBUG
        std::cout << "\tCC dfs:\n";
#endif
        while (dfs.size()) {
            auto v = dfs.back();
            dfs.pop_back();
            if (visited.count(v))   continue;
            visited.insert(v);
#ifdef DEBUG
            auto vd = v->detector;
            std::cout << "\t\t" << vd.first << "(" << c2i(vd.second) << ") to";
#endif
            for (auto w : connection_graph.get_neighbors(v)) {
                if (w == prev[v])   continue;
                cdet_t wd = w->detector;
#ifdef DEBUG
                std::cout << " " << wd.first << "(" << c2i(wd.second) << ")";
#endif
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
                    __COLOR cc_color;
                    if (wd.second == __COLOR::red || wd.second == __COLOR::blue)  {
                        cc_color = __COLOR::green;
                    } else {
                        cc_color = __COLOR::blue;
                    }
                    jamming_set_t js = std::make_pair(cc, cc_color);
                    jamming_sets.push_back(js);
                } else if (!visited.count(w)) {
                    prev[w] = v;
                    dfs.push_back(w);
                }
            }
#ifdef DEBUG
            std::cout << "\n";
#endif
        }
    }
    // Add everything not in a connected components to its own jamming set.
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
    jamming_sets.push_back(std::make_pair(not_cc, __COLOR::red));
    // Compute correction using the jamming sets.
    stim::simd_bits corr(n_observables);
    corr.clear();
#ifdef DEBUG
    std::cout << "Jamming sets:\n";
#endif
    for (auto& p1 : jamming_sets) {
        // We track the incident colors for each detector
        std::map<cdet_t, __COLOR> detector_to_incident_color;
        auto& js = p1.first;
        __COLOR jsc = p1.second;
#ifdef DEBUG
        std::cout << "\tColor: " << c2i(jsc) << "\n";
#endif
        for (auto p2 : js) {
            cdet_t cd1 = std::get<0>(p2);
            cdet_t cd2 = std::get<1>(p2);
            __COLOR restricted_color = std::get<2>(p2);
#ifdef DEBUG
            std::cout << "\t\t" << cd1.first << "(" << c2i(cd1.second) << ") <----> "
                << cd2.first << "(" << c2i(cd2.second) << ")\t COLOR = " 
                << c2i(restricted_color) << "\n";
#endif
            if (restricted_color == jsc) {
                continue;
            }
            // Now get decoding graph.
            const uint dec_index = 2 - c2i(restricted_color);
            DecodingGraph& gr = rlatt_dec[dec_index]->decoding_graph;
            uint rlatt_d1 = cd1.first == BOUNDARY_INDEX ? 
                                BOUNDARY_INDEX : to_rlatt[cd1.first][dec_index];
            uint rlatt_d2 = cd2.first == BOUNDARY_INDEX ? 
                                BOUNDARY_INDEX : to_rlatt[cd2.first][dec_index];
            auto v1 = gr.get_vertex(rlatt_d1);
            auto v2 = gr.get_vertex(rlatt_d2);
            // Perform frame changes along the path provided one of the detectors
            // has the same color as the jamming set.
#ifdef DEBUG
            std::cout << "\t\t\tpath (length = "
                << gr.get_error_chain_data(v1, v2).chain_length<< "):\n";
#endif
            auto path = gr.get_error_chain_data(v1, v2).error_chain;
            __COLOR possible_boundary_color1, possible_boundary_color2;
            if (restricted_color == __COLOR::red) {
                possible_boundary_color1 = __COLOR::blue;
                possible_boundary_color2 = __COLOR::green;
            } else if (restricted_color == __COLOR::blue) {
                possible_boundary_color1 = __COLOR::red;
                possible_boundary_color2 = __COLOR::green;
            } else {
                possible_boundary_color1 = __COLOR::red;
                possible_boundary_color2 = __COLOR::blue;
            }
            for (uint i = 1; i < path.size(); i++) {
                auto vx = path[i-1];
                auto vy = path[i];
                
                uint64_t dx = vx->id;
                uint64_t dy = vy->id;
                if (dx == BOUNDARY_INDEX)   std::swap(dx, dy);
                dx = from_rlatt[std::make_pair(dx, dec_index)];
                __COLOR dx_color = color_map[dx];
                __COLOR dy_color;
                if (dy == BOUNDARY_INDEX) {
                    if (dx_color == possible_boundary_color1) {
                        dy_color = possible_boundary_color2;
                    } else {
                        dy_color = possible_boundary_color1;
                    }
                } else {
                    dy = from_rlatt[std::make_pair(dy, dec_index)];
                    dy_color = color_map[dy];
                }
#ifdef DEBUG
                std::cout << "\t\t\t\tE : " << dx << "(" << c2i(dx_color) << ") , " 
                    << dy << "(" << c2i(dy_color) << ")";
#endif
                if (dx_color != jsc && dy_color != jsc) {
#ifdef DEBUG
                    std::cout << "S\n";
#endif
                    continue;
                }
                auto e = gr.get_edge(vx, vy);
                if (e->frames.empty() && dx != BOUNDARY_INDEX && dy != BOUNDARY_INDEX) {
#ifdef DEBUG
                    std::cout << "E\n";
#endif
                    continue;
                }
                cdet_t cdx = std::make_pair(dx, dx_color);
                cdet_t cdy = std::make_pair(dy, dy_color);
                __COLOR old_dx_incident_color = detector_to_incident_color.count(cdx)
                                                    ? detector_to_incident_color[cdx]
                                                    : __COLOR::none;
                __COLOR old_dy_incident_color = detector_to_incident_color.count(cdy)
                                                    ? detector_to_incident_color[cdy]
                                                    : __COLOR::none;
                detector_to_incident_color[cdx] = restricted_color;
                detector_to_incident_color[cdy] = restricted_color;
                if (old_dx_incident_color != __COLOR::none) {
                    if (old_dx_incident_color != restricted_color) {
#ifdef DEBUG
                        std::cout << "\tJ\n";
#endif
                        continue;
                    }
                }
                if (old_dy_incident_color != __COLOR::none) {
                    if (old_dy_incident_color != restricted_color) {
#ifdef DEBUG
                        std::cout << "\tJ\n";
#endif
                        continue;
                    }
                }
                for (auto f : e->frames) {
                    corr[f] ^= 1;                
#ifdef DEBUG
                    std::cout << "\tFR" << f;
#endif
                }
#ifdef DEBUG
                std::cout << "\n";
#endif
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
#ifdef DEBUG
    std::cout << "\tDetectors on lattice restricting color "
            << c2i(restricted_color) << ":";
#endif
    for (uint d : detectors) {
        if (color_map[d] == restricted_color)   continue;
        uint dd = to_rlatt[d][i];
        restricted_syndrome[dd] = 1;
#ifdef DEBUG
        std::cout << " " << d;
#endif
    }
#ifdef DEBUG
    std::cout << "\n";
#endif
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

}   // decoder
}   // qontra
