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

#define CDET_TO_ID(x)   ((x.first) | ((uint64_t)c2i(x.second)) << 60)

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
        uint    rlatt;
    };
    typedef Graph<cv_t, ce_t>   ConnGraph;

    ConnGraph connection_graph;
    Decoder::result_t* res_array[] = { &r_lattRG, &r_lattRB, &r_lattGB };

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
            uint64_t id1 = CDET_TO_ID(cd1);
            uint64_t id2 = CDET_TO_ID(cd2);
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
            e->rlatt = i;
            connection_graph.add_edge(e);
        }
    }
    // Here, we will partition the syndrome into two parts:
    //  (1) Syndrome bits not in a connected component with the red boundary.
    //  (2) Everything else.
    // Then, we will "jam" matchings in each partition separately. Jamming is essentially
    // when a detector has matchings in two different colors (i.e. G and B). Then, we ignore
    // one of the incident frame changes as this implies that we are at a face.
    std::set<cdet_t> in_connected_components;

    typedef std::pair<std::vector<cdet_t>, __COLOR> jamming_set_t;
    
    cdet_t red_boundary = std::make_pair(BOUNDARY_INDEX, __COLOR::red);
    auto rbv = connection_graph.get_vertex(CDET_TO_ID(red_boundary));

    std::vector<jamming_set_t> jamming_sets;
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
            visited.insert(v);
            for (auto w : connection_graph.get_neighbors(v)) {
                cdet_t wd = w->detector;
                if (wd.first == BOUNDARY_INDEX) {
                    // This is one connected component.
                    std::vector<cdet_t> cc;
                    auto curr = v;
                    while (curr != rbv) {
                        cc.push_back(v->detector);
                        curr = prev[curr];
                    }
                    cc.push_back(wd);
                    cc.push_back(red_boundary);
                    __COLOR cc_color;
                    if (wd.second == __COLOR::red || wd.second == __COLOR::blue)  {
                        cc_color = __COLOR::green;
                    } else {
                        cc_color = __COLOR::blue;
                    }
                    jamming_set_t js = std::make_pair(cc, cc_color);
                    for (auto x : cc)   in_connected_components.insert(x);
                    jamming_sets.push_back(js);
                } else if (!visited.count(w)) {
                    prev[w] = v;
                    dfs.push_back(w);
                }
            }
        }
    }
    // Add everything not in a connected components to its own jamming set.
    std::vector<cdet_t> not_cc;
    for (uint d : detectors) {
        __COLOR c = color_map[d];
        cdet_t det = std::make_pair(d, c);
        if (!in_connected_components.count(det))    not_cc.push_back(det);
    }
    jamming_sets.push_back(std::make_pair(not_cc, __COLOR::red));
    // Compute correction using the jamming sets.
    stim::simd_bits corr(n_observables);
    corr.clear();

#ifdef DEBUG
    std::cout << "Jamming sets:\n";
#endif
    for (auto& pair : jamming_sets) {
        std::map<uint, __COLOR> detector_to_color;
        auto& js = pair.first;
        __COLOR jsc = pair.second;
#ifdef DEBUG
        std::cout << "\tColor: " << c2i(jsc) << "\n";
#endif
        for (uint i = 1; i < js.size(); i++) {
            cdet_t cd1 = js[i-1];
            cdet_t cd2 = js[i];
#ifdef DEBUG
            std::cout << "\t\t" << cd1.first << "(" << c2i(cd1.second) << ") <----> "
                << cd2.first << "(" << c2i(cd2.second) << ")\n";
#endif
            // Get the path between cdx and cdy and only examine detectors
            // with color jsc.
            //
            // First, check the colors of cdx and cdy and then get the corresponding
            // decoding graph.
            __COLOR restricted_color;
            // Compute restricted color.
            std::set<__COLOR> available{__COLOR::red, __COLOR::green, __COLOR::blue};
            __COLOR c1 = cd1.second;
            __COLOR c2 = cd2.second;
            available.erase(c1);
            available.erase(c2);
            restricted_color = *available.begin();
            // Now get decoding graph.
#ifdef DEBUG
            std::cout << "\t\t\trestricted = " << c2i(restricted_color) << "\n";
#endif
            const uint dec_index = 2 - c2i(restricted_color);
            DecodingGraph& gr = rlatt_dec[dec_index]->decoding_graph;
            uint rlatt_d1 = to_rlatt[cd1.first][dec_index];
            uint rlatt_d2 = to_rlatt[cd2.first][dec_index];
            auto v1 = gr.get_vertex(rlatt_d1);
            auto v2 = gr.get_vertex(rlatt_d2);
            // Get path between v1 and v2: only perform frame changes incident
            // on vertices matching the color jsc.
            auto error_chain = gr.get_error_chain_data(v1, v2).error_chain;
            for (uint j = 1; j < error_chain.size(); j++) {
                auto vx = error_chain[j-1];
                auto vy = error_chain[j];
                uint rlatt_dx = vx->id;
                uint rlatt_dy = vy->id;
                uint dx = from_rlatt[std::make_pair(rlatt_dx, dec_index)];
                uint dy = from_rlatt[std::make_pair(rlatt_dy, dec_index)];
#ifdef DEBUG
                std::cout << "\t\t\t(" << dx << ", " << dy << ")(";
#endif
                if (color_map[dx] == jsc) {
                    auto old_color = detector_to_color.count(dx) 
                                        ? detector_to_color[dx] : restricted_color;
                    detector_to_color[dx] = restricted_color;
                    // Jam the correction if necessary.
                    if (old_color != restricted_color) {
#ifdef DEBUG
                        std::cout << "j)\n";
#endif
                        continue;
                    }
                } else {
                    auto old_color = detector_to_color.count(dy)
                                        ? detector_to_color[dy] : restricted_color;
                    detector_to_color[dy] = restricted_color;
                    // Jam the correction if necessary.
                    if (old_color != restricted_color) {
#ifdef DEBUG
                        std::cout << "j)\n";
#endif
                        continue;
                    }
                }
                auto e = gr.get_edge(vx, vy);
                for (auto fr : e->frames) {
                    corr[fr] ^= 1;
#ifdef DEBUG
                    std::cout << fr;
#endif
                }
#ifdef DEBUG
                std::cout << ")\n";
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

    stim::simd_bits restricted_syndrome(n_detectors+n_observables);
    restricted_syndrome.clear();
    uint hw = 0;
    for (uint d : detectors) {
        if (color_map[d] == restricted_color)   continue;
        uint dd = to_rlatt[d][i];
        restricted_syndrome[d] = 1;
        hw++;
    }
    auto res = dec->decode_error(restricted_syndrome);
    // Convert the error assignments to the original circuit's detectors.
    for (auto& aa : res.error_assignments) {
        uint& d1 = std::get<0>(aa);
        uint& d2 = std::get<1>(aa);
        if (d1 != BOUNDARY_INDEX)   d1 = from_rlatt[std::make_pair(d1, i)];
        if (d2 != BOUNDARY_INDEX)   d2 = from_rlatt[std::make_pair(d2, i)];
    }
    return res;
}

}   // decoder
}   // qontra
