/*
 *  author: Suhas Vittal
 *  date:   19 July 2023
 * */

#include "decoder/restriction.h"

#include <assert.h>

#define DEBUG

namespace qontra {
namespace decoder {

namespace restriction {

void
inherit_neighbors(StructureGraph& gr, stv_t* v, stv_t* w, std::function<bool(stv_t*)> inherit_cb) {
    for (auto u : gr.get_neighbors(w)) {
        if (u == v || gr.contains(v, u))    continue;
        ste_t* e = new ste_t;
        e->src = v;
        e->dst = u;
        e->is_undirected = true;
        gr.add_edge(e);
        if (inherit_cb(u))   inherit_neighbors(gr, v, u, inherit_cb);
    }
}

}   // restriction



using namespace graph;
using namespace restriction;

// TODO when finally adding measurement errors, do the following:
//  (1) Modify the constructor to request detectors per round.
//  (2) Modify StructureGraph construction to avoid creating new vertices
//      each round.
//  (3) Modify get_all_edges_in_component to flatten measurement results.

RestrictionDecoder::RestrictionDecoder(const stim::Circuit& circuit)
    :Decoder(circuit, false),
    lattice_structure(),
    boundary_adjacent(),
    rlatt_dec(),
    to_rlatt(),
    from_rlatt(),
    color_map()
{
    boundary_adjacent.fill(std::set<cdet_t>());
    // Restricted lattice order: RG (so blue is restricted), RB, GB
    const __COLOR colors[3] = { __COLOR::blue, __COLOR::green, __COLOR::red };
    // Create new Stim circuits for each restricted color.
    //
    // Also build the StructureGraph. We do so by infering connectivity
    // from the CNOT operations. Initially the IDs will be the qubit numbers
    // in the Stim circuit, but we will convert them to match the colored
    // detector values.
    std::set<stv_t*> unused; // We will delete all vertices in this set afterwards.
    for (uint i = 0; i < circuit.count_qubits(); i++) {
        stv_t* v = new stv_t;
        v->id = i;
        v->detector = std::make_pair(i, __COLOR::none);
        unused.insert(v);
        lattice_structure.add_vertex(v);
    }
    std::deque<stv_t*> measurement_queue;

    uint det_ctr = 0;

    std::array<uint, 3> subdet_ctr;
    subdet_ctr.fill(0);

    std::array<stim::Circuit, 3> subcircuits;
    subcircuits.fill(stim::Circuit());

    circuit.for_each_operation([&] (stim::Operation op) {
        std::string gate_name(op.gate->name);
        if (gate_name == "DETECTOR") {
            int color_id = (int) op.target_data.args[0];
            // Add to subcircuit
            for (uint i = 0; i < 3; i++) {
                __COLOR restricted_color = colors[i];
                stim::Circuit& subcircuit = subcircuits[i];
                if (i2c(color_id) != restricted_color) {
                    subcircuit.append_operation(op);

                    to_rlatt[det_ctr][i] = subdet_ctr[i];
                    from_rlatt[std::make_pair(subdet_ctr[i], i)] = det_ctr;
                    subdet_ctr[i]++;
                }
            }
            // Modify vertex in structure graph and remove the vertex from the
            // unused set.
            uint64_t offset = op.target_data.targets[0].data & ~stim::TARGET_RECORD_BIT;
            stv_t* v = measurement_queue[measurement_queue.size() - offset];
            cdet_t cd = std::make_pair(det_ctr, i2c(color_id));
            v->detector = cd;
            lattice_structure.change_id(v, CDET_TO_ID(cd));
            unused.erase(v);
            color_map[det_ctr++] = i2c(color_id);
            return;
        } else if (gate_name == "M" || gate_name == "MR") {
            const auto& targets = op.target_data.targets;
            for (uint i = 0; i < targets.size(); i++) {
                auto v = lattice_structure.get_vertex(targets[i].data);
                measurement_queue.push_back(v); 
            }
        } else if (gate_name == "CX") {
            const auto& targets = op.target_data.targets;
            for (uint i = 0; i < targets.size(); i += 2) {
                uint q1 = targets[i].data;
                uint q2 = targets[i+1].data;
                auto v1 = lattice_structure.get_vertex(q1);
                auto v2 = lattice_structure.get_vertex(q2);
                if (!lattice_structure.contains(v1, v2)) {
                    ste_t* e = new ste_t;
                    e->src = v1;
                    e->dst = v2;
                    e->is_undirected = true;
                    lattice_structure.add_edge(e);
                }
            }
        }
        for (uint i = 0; i < 3; i++) {
            subcircuits[i].append_operation(op);
        }
    });
    // If two used vertices are already connected, then one of them is a gauge and 
    // both inherit each others unused vertices: we must check for this first. 
    //
    // Now, we need to create edges between the used vertices in the StructureGraph.
    // We do this by connecting two vertices if they share another vertex in common.
    //
    // We will also add three boundaries and connect used vertices
    // if necessary.
    for (auto v : lattice_structure.get_vertices()) {
        if (unused.count(v))    continue;
        for (auto w : lattice_structure.get_neighbors(v)) {
            if (unused.count(w))    continue;
            inherit_neighbors(lattice_structure, v, w, [&] (stv_t* u) {
                return !unused.count(u);
            });
        }
    }

    // Add boundaries and connect them to one another.
    for (uint i = 0; i < 3; i++) {
        __COLOR boundary_color = i2c(i);
        cdet_t cd = std::make_pair(BOUNDARY_INDEX, boundary_color);
        stv_t* v = new stv_t;
        v->id = CDET_TO_ID(cd);
        v->detector = cd;
        lattice_structure.add_vertex(v);
    }
    for (uint i = 0; i < 3; i++) {
        cdet_t x = std::make_pair(BOUNDARY_INDEX, i2c(i));
        auto bv = lattice_structure.get_vertex(CDET_TO_ID(x));
        for (uint j = i+1; j < 3; j++) {
            cdet_t y = std::make_pair(BOUNDARY_INDEX, i2c(i));
            auto bw = lattice_structure.get_vertex(CDET_TO_ID(y));
            ste_t* e = new ste_t;
            e->src = bv;
            e->dst = bw;
            e->is_undirected = true;
            lattice_structure.add_edge(e);

            boundary_adjacent[i].insert(y);
            boundary_adjacent[j].insert(x);
        }
    }

    for (auto v : lattice_structure.get_vertices()) {
        if (unused.count(v))    continue;
        if (v->detector.first == BOUNDARY_INDEX)    continue;
        auto tmp = lattice_structure.get_neighbors(v);
        std::set<stv_t*> vadj(tmp.begin(), tmp.end());
        for (auto w : lattice_structure.get_vertices()) {
            if (unused.count(w))    continue;
            if (v == w || lattice_structure.contains(v, w)) continue;
            if (w->detector.first == BOUNDARY_INDEX)    continue;
            tmp = lattice_structure.get_neighbors(w);
            std::set<stv_t*> wadj(tmp.begin(), tmp.end());
            
            bool intersects = false;
            for (auto x : vadj) {
                if (wadj.count(x) && unused.count(x)) {
                    intersects = true;
                    break;
                }
            }
            if (intersects) {
                ste_t* e = new ste_t;
                e->src = v;
                e->dst = w;
                e->is_undirected = true;
                if (!lattice_structure.add_edge(e)) std::cout << "\t\tHUHHH?\n";
            }
        }
        
        // Check if we need to add a boundary connection. This is only required if
        // the number of adjacent detectors is less than the number of adjacent data
        // qubits (that is, adjacent used < adjacent unused).
        uint number_of_adjacent_unused = 0;
        std::array<uint, 3> number_of_adjacent_used;    // one for each color.
        number_of_adjacent_used.fill(0);
        for (auto w : lattice_structure.get_neighbors(v)) {
            if (unused.count(w))    number_of_adjacent_unused++;
            else {
                int i = c2i(w->detector.second);
                number_of_adjacent_used[i]++;
            }
        }
        uint sum_number_of_adjacent_used = 
            number_of_adjacent_used[0] + number_of_adjacent_used[1] + number_of_adjacent_used[2];
        if (number_of_adjacent_unused > sum_number_of_adjacent_used) {
            __COLOR zero_color = v->detector.second;    // As expected, we should not be
                                                        // adjacent to any detectors of the same
                                                        // color.
            const uint half_sum = number_of_adjacent_unused >> 1;
            for (uint i = 0; i < 3; i++) {
                __COLOR c = i2c(i);
                if (c == zero_color)    continue;
                if (number_of_adjacent_used[i] < half_sum) {
                    cdet_t cb = std::make_pair(BOUNDARY_INDEX, c);
                    auto wb = lattice_structure.get_vertex(CDET_TO_ID(cb));
                    ste_t* e = new ste_t;
                    e->src = v;
                    e->dst = wb;
                    e->is_undirected = true;
                    lattice_structure.add_edge(e);
                    boundary_adjacent[i].insert(v->detector);
                }
            }
        }
    }
    // Delete all unused.
    for (auto x : unused) {
        lattice_structure.delete_vertex(x);
    }
#ifdef DEBUG
    std::cout << "***Lattice structure:\n";
    for (auto v : lattice_structure.get_vertices()) {
        auto cd1 = v->detector;
        std::cout << cd1.first << "(" << c2i(cd1.second) << "):";
        for (auto w : lattice_structure.get_neighbors(v)) {
            auto cd2 = w->detector;
            std::cout << " " << cd2.first << "(" << c2i(cd2.second) << ")";
        }
        std::cout << "\n";
    }
#endif

    for (uint i = 0; i < 3; i++) {
        rlatt_dec[i] = new MWPMDecoder(subcircuits[i]);
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
            e->color = restricted_colors[i];
            connection_graph.add_edge(e);
        }
    }
    // Here, we will partition the syndrome into two parts:
    //  (1) Syndrome bits not in a connected component with the red boundary.
    //  (2) Everything else.
    cdet_t red_boundary = std::make_pair(BOUNDARY_INDEX, __COLOR::red);
    auto rbv = connection_graph.get_vertex(CDET_TO_ID(red_boundary));

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

        auto v = connection_graph.get_vertex(CDET_TO_ID(det));
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
                cdet_t cdz;
                cdetpair_t e_xy, e_xz, e_yz, e_yx, e_zx, e_zy;
                for (auto it = common.begin(); it != common.end(); it++) {
                    cdz = *it;
                    e_xy = std::make_pair(cdx, cdy);
                    e_xz = std::make_pair(cdx, cdz);
                    e_yz = std::make_pair(cdy, cdz);

                    e_yx = std::make_pair(cdy, cdx);
                    e_zx = std::make_pair(cdz, cdx);
                    e_zy = std::make_pair(cdz, cdy);

                    if (!(in_face.count(e_xy) || in_face.count(e_yx)
                        || in_face.count(e_xz) || in_face.count(e_zx)
                        || in_face.count(e_yz) || in_face.count(e_zy)))
                    {
                        goto cdz_is_valid;
                    }
                }
                break;
cdz_is_valid:
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
                // Check if d1 is even adjacent to this boundary.
                cdet_t cb1 = std::make_pair(d2, c2);
                if (!boundary_adjacent[c2i(c2)].count(cd1)) {
                    // Then, this path must pass through another boundary.
                    int option1 = (c2i(c2) + 1) % 3;
                    int option2 = (c2i(c2) + 2) % 3;
                    __COLOR b2c;
                    if (boundary_adjacent[option1].count(cd1)) {
                        b2c = i2c(option1);
                    } else {
                        b2c = i2c(option2);
                    }
                    cdet_t cb2 = std::make_pair(BOUNDARY_INDEX, b2c);
                    // Get common neighbor between the two boundaries.
                    auto common = get_common_neighbors(cb1, cb2);
                    cdet_t idet = *common.begin();
                    edges.insert(std::make_pair(cd1, cb2));
                    edges.insert(std::make_pair(cb2, idet));
                    edges.insert(std::make_pair(idet, cb1));
                    continue;
                }
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
    stv_t* v = lattice_structure.get_vertex(CDET_TO_ID(cd));
    for (auto w : lattice_structure.get_neighbors(v)) {
        neighbors.insert(w->detector);
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
