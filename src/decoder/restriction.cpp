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
using namespace restriction;

// TODO when finally adding measurement errors, do the following:
//  (1) Modify the constructor to request detectors per round.
//  (2) Modify StructureGraph construction to avoid creating new vertices
//      each round.
//  (3) Modify get_all_edges_in_component to flatten measurement results.

RestrictionDecoder::RestrictionDecoder(
        const stim::Circuit& circuit,
        const uint detectors_per_round,
        StructureGraph* stgr, 
        const std::map<uint64_t, uint64_t>& d2b,
        const std::map<uint64_t, uint64_t>& f2o)
    :Decoder(circuit, false),
    lattice_structure(stgr),
    lattice_matrix(),
    boundary_adjacent(),
    rlatt_dec(),
    to_rlatt(),
    from_rlatt(),
    color_map(),
    detectors_per_round(detectors_per_round),
    detector_to_base(d2b),
    flag_to_owner(f2o)
{
    // Build the lattice structure matrix.
    ewf_t<stv_t> d_w = [&] (stv_t* x, stv_t* y) {
        return (x->is_flag || y->is_flag) ? 10000000 : 1;
    };
    distance::callback_t<stv_t, std::vector<stv_t*>> d_cb = 
        [] (stv_t* v1, 
                stv_t* v2, 
                const std::map<stv_t*, fp_t>& dist,
                const std::map<stv_t*, stv_t*>& pred)
        {
            std::vector<stv_t*> path;
            auto curr = v2;
            while (curr != v1) {
                path.push_back(curr);
                curr = pred.at(curr);
            }
            path.push_back(v1);
            std::reverse(path.begin(), path.end());
            return path;
        };
    lattice_matrix = distance::create_distance_matrix(stgr, d_w, d_cb);
    detector_to_base[BOUNDARY_INDEX] = BOUNDARY_INDEX;
    boundary_adjacent.fill(std::set<cdet_t>());
    // Restricted lattice order: RG (so blue is restricted), RB, GB
    const __COLOR colors[3] = { __COLOR::blue, __COLOR::green, __COLOR::red };
    // Create new Stim circuits for each restricted color.
    uint det_ctr = 0;
    std::array<uint, 3> subdet_ctr;
    std::array<stim::Circuit, 3> subcircuits;

    subdet_ctr.fill(0);
    subcircuits.fill(stim::Circuit());

    circuit.for_each_operation([&] (stim::Operation op) {
        std::string gate_name(op.gate->name);
        if (gate_name == "DETECTOR") {
            int color_id = (int) op.target_data.args[0];
            // Add to subcircuit
            for (uint dec_index = 0; dec_index < 3; dec_index++) {
                __COLOR restricted_color = colors[dec_index];
                stim::Circuit& subcircuit = subcircuits[dec_index];
                if (i2c(color_id) != restricted_color) {
                    subcircuit.append_operation(op);

                    to_rlatt[det_ctr][dec_index] = subdet_ctr[dec_index];
                    from_rlatt[std::make_pair(subdet_ctr[dec_index], dec_index)] = det_ctr;
                    subdet_ctr[dec_index]++;
                }
            }
            color_map[det_ctr] = i2c(color_id);
            det_ctr++;
        } else {
            for (uint dec_index = 0; dec_index < 3; dec_index++) {
                subcircuits[dec_index].append_operation(op);
            }
        }
    });

    for (uint dec_index = 0; dec_index < 3; dec_index++) {
        rlatt_dec[dec_index] = new MWPMDecoder(subcircuits[dec_index]);
    }

#ifdef DEBUG
    std::cout << "***Lattice structure:\n";
    for (auto v : lattice_structure->get_vertices()) {
        auto cd1 = v->detector;
        std::cout << cd1.first << "(" << c2i(cd1.second) << "):";
        for (auto w : lattice_structure->get_neighbors(v)) {
            auto cd2 = w->detector;
            std::cout << " " << cd2.first << "(" << c2i(cd2.second) << ")";
        }
        std::cout << "\n";
    }
#endif
    // Check boundary adjacency.
    for (uint dec_index = 0; dec_index < 3; dec_index++) {
        cdet_t boundary = std::make_pair(BOUNDARY_INDEX, i2c(dec_index));
        auto bv = lattice_structure->get_vertex(CDET_TO_ID(boundary));
        for (auto bw : lattice_structure->get_vertices()) {
            if (lattice_structure->contains(bv, bw)) {
                boundary_adjacent[dec_index].insert(bw->detector);
            }
        }
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
#endif
    for (auto x : detectors) {
#ifdef DEBUG
        std::cout << " " << x << "(" << c2i(color_map[x]) << ")";
#endif
    }
#ifdef DEBUG
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
    for (uint dec_index = 0; dec_index < 3; dec_index++) {
        auto res_p = res_array[dec_index];
        for (auto aa : res_p->error_assignments) {
            uint d1 = std::get<0>(aa);
            uint d2 = std::get<1>(aa);

            if (d1 == BOUNDARY_INDEX)   std::swap(d1, d2);
            __COLOR c1 = color_map[d1];
            // Compute d2's color. If d2 is the boundary, then we need
            // to use d1's color to do so.
            __COLOR c2;
            if (d2 == BOUNDARY_INDEX) {
                uint rlatt_d1 = to_rlatt[d1][dec_index];
                // Get decoding graph vertices.
                DecodingGraph& gr = rlatt_dec[dec_index]->decoding_graph;
                auto vd1 = gr.get_vertex(rlatt_d1);
                auto vd2 = gr.get_vertex(BOUNDARY_INDEX);
                auto ec_data = gr.get_error_chain_data(vd1, vd2);
                // If the error chain length is even, then d2 is the same
                // color as d1. Otherwise, it is the opposite color.
                uint len = ec_data.chain_length;
                if (len & 0x1) {
                    c2 = get_remaining_color(c1, restricted_colors[dec_index]);
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
            e->color = restricted_colors[dec_index];
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
    connected_components.push_back(std::make_pair(not_cc, __COLOR::none));
                                        // Label it as none so we know it is "not cc".
    stim::simd_bits corr(n_observables);
    corr.clear();

#ifdef DEBUG
    std::cout << "Connected Components:\n";
#endif
    std::set<cdetpair_t> in_face_cc;
    std::set<cdetpair_t> in_face_not_cc;

    for (auto& pair : connected_components) {
        auto cc = pair.first;
        auto& cc_color = pair.second;

        std::set<cdetpair_t>& in_face = cc_color == __COLOR::none ? 
                                            in_face_not_cc : in_face_cc;
#ifdef DEBUG
        std::cout << "\tColor " << c2i(cc_color) << ":\n";
#endif
        if (cc_color == __COLOR::none) {
            cc_color = __COLOR::red;
        }
        
        auto flipped_edges = get_all_edges_in_component(pair);
        std::set<cdet_t> vertices_in_flips; // of the same color as the component.
        std::map<cdet_t, std::set<cdet_t>> incident_vertices;
#ifdef DEBUG
        std::cout << "\t\tEdges:\n";
#endif
        for (cdetpair_t edge : flipped_edges) { 
            cdet_t cdx = edge.first;
            cdet_t cdy = edge.second;
#ifdef DEBUG
            std::cout << "\t\t\t" << cdx.first << "(" << c2i(cdx.second) << ")"
                        << ", " << cdy.first << "(" << c2i(cdy.second) << ")\n";
#endif
            if (cdx.second == cdy.second)   continue;
            if (cdx.second == cc_color) {
                vertices_in_flips.insert(cdx);
                incident_vertices[cdx].insert(cdy);
            } 
            if (cdy.second == cc_color) {
                vertices_in_flips.insert(cdy);
                incident_vertices[cdy].insert(cdx);
            }
        }
        // Now, we will apply corrections on faces around each vertex.
        for (cdet_t cdx : vertices_in_flips) {
#ifdef DEBUG
            std::cout << "\t\tChecking vertex " << cdx.first << "(" << c2i(cdx.second) << "):\n";
#endif
            // Get all the faces incident to cdx.
            typedef std::tuple<cdet_t, cdet_t, cdet_t>  face_t;
            std::set<face_t> all_faces;
            auto neighbors = get_neighbors(cdx);
#ifdef DEBUG
            std::cout << "\t\t\tneighbors:";
            for (auto x : neighbors) std::cout << " " << x.first << "(" << c2i(x.second) << ")";
            std::cout << "\n";
#endif
        
#ifdef DEBUG
            std::cout << "\t\t\tfaces:\n";
#endif
            for (auto cdy : neighbors) {
                // Get common neighbor of cdx and cdy.
                auto n = get_neighbors(cdy);
                std::set<cdet_t> common;
                std::set_intersection(neighbors.begin(), neighbors.end(),
                                        n.begin(), n.end(),
                                        std::inserter(common, common.begin()));
                for (auto it = common.begin(); it != common.end(); it++) {
                    cdet_t cdz = *it;
                    if (cdz == cdy) continue;
                    face_t f0 = std::make_tuple(cdx, cdy, cdz);
                    face_t f1 = std::make_tuple(cdx, cdz, cdy); // Make sure we are not double
                                                                // counting.
                    if (!all_faces.count(f0) && !all_faces.count(f1)) {
                        all_faces.insert(f0);
#ifdef DEBUG
                        std::cout << "\t\t\t\t< " << cdx.first << "(" << c2i(cdx.second) << "), "
                                    << cdy.first << "(" << c2i(cdy.second) << "), "
                                    << cdz.first << "(" << c2i(cdz.second) << ")\n";
#endif
                    }
                }
            }
            // Now, compute the best subset of faces via exhaustive search.
            uint64_t max_satisfied_edges = 0;
            uint64_t best_ctr_num_faces = 0;
            uint64_t best_ctr = 0;
            uint64_t ctr = 0;
            const uint64_t max_ctr = 1 << all_faces.size();
            while (ctr < max_ctr) {
                std::set<cdetpair_t> boundary_edges;    
                        // The edges on the boundary of the face subset.
                uint64_t dec_index = 0;
                uint64_t num_faces_used = 0;
                for (auto it = all_faces.begin(); it != all_faces.end(); it++) {
                    if (!(ctr & (1 << (dec_index++)))) {
                        continue;
                    }
                    num_faces_used++;
                    face_t f = *it;
                    cdet_t x = std::get<0>(f), y = std::get<1>(f), z = std::get<2>(f);
                    auto xy = std::make_pair(x, y);
                    auto xz = std::make_pair(x, z);
                    auto yz = std::make_pair(y, z);
                    std::vector<cdetpair_t> edges_of_face{xy, xz, yz};
                    for (auto e : edges_of_face) {
                        xor_pair_into(boundary_edges, e);
                    }
                }
                // Now, check how many flipped edges are on the boundary.
                // If any edges are already in a face, then this is invalid.
                uint64_t satisfied_edges = 0;
                bool any_in_face = false;
                bool any_incident_not_flipped = false;
                for (auto e : boundary_edges) {
                    any_in_face |= in_face.count(e);
                    if (e.first == cdx || e.second == cdx) {
                        any_incident_not_flipped |= !flipped_edges.count(e)
                                                        && !(e.first.first == BOUNDARY_INDEX
                                                                && e.second.first == BOUNDARY_INDEX);
                    }
                    satisfied_edges += flipped_edges.count(e);
                }
                if (any_in_face || any_incident_not_flipped) {
                    ctr++;
                    continue;
                }
                if (satisfied_edges > max_satisfied_edges
                    || (satisfied_edges == max_satisfied_edges
                        && num_faces_used < best_ctr_num_faces)) 
                {
                    best_ctr = ctr;
                    max_satisfied_edges = satisfied_edges;
                    best_ctr_num_faces = num_faces_used;
                }
                ctr++;
            }
            // Now, update the correction as well as which edges are in a face.
            if (best_ctr == 0) {
#ifdef DEBUG
                std::cout << "\t\t\tno good local correction.\n";
#endif
                continue;
            }
            uint dec_index = 0;
            for (auto it = all_faces.begin(); it != all_faces.end(); it++) {
                if (!(best_ctr & (1 << (dec_index++)))) continue;
                face_t f = *it;
                cdet_t x = std::get<0>(f), y = std::get<1>(f), z = std::get<2>(f);
                corr ^= get_correction_for_face(x, y, z);
                auto xy = std::make_pair(x, y);
                auto yx = std::make_pair(y, x);
                auto xz = std::make_pair(x, z);
                auto zx = std::make_pair(z, x);
                auto yz = std::make_pair(y, z);
                auto zy = std::make_pair(z, y);
                std::vector<cdetpair_t> edges_of_face{xy, yx, xz, zx, yz, zy};
                for (auto e : edges_of_face) {
                    if (flipped_edges.count(e)) insert_pair_into(in_face, e);
                }
            }
        }
    }

#ifdef DEBUG
    std::cout << "\tobs =";
    for (uint dec_index = 0; dec_index < n_observables; dec_index++) {
        std::cout << syndrome[n_detectors+dec_index]+0;
        if (dec_index == 0)  std::cout << "|";
    }
    std::cout << "\tcorr =";
    for (uint dec_index = 0; dec_index < n_observables; dec_index++) {
        std::cout << corr[dec_index]+0;
        if (dec_index == 0)  std::cout << "|";
    }
    std::cout << "\n";
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
    const uint dec_index = 2 - c2i(restricted_color);
    MWPMDecoder* dec = rlatt_dec[dec_index];
    DecodingGraph& gr = dec->decoding_graph;

    const uint n_detectors = dec->get_circuit().count_detectors();
    const uint n_observables = dec->get_circuit().count_observables();

    syndrome_t restricted_syndrome(n_detectors+n_observables);
    restricted_syndrome.clear();

    std::set<std::pair<decoding::vertex_t*, decoding::vertex_t*>> updated_edges;
    std::vector<decoding::vertex_t*> vertices_in_syndrome;
    for (uint di : detectors) {
        if (color_map[di] == restricted_color)   continue;
        uint ddi = to_rlatt[di][dec_index];
        // Check if d is a flag. If so, make changes to the weights of the decoder.
        cdet_t cdi = std::make_pair(detector_to_base[di], color_map[di]); 
        stv_t* lattice_v = lattice_structure->get_vertex(CDET_TO_ID(cdi));

        decoding::vertex_t* dv = gr.get_vertex(ddi);
        if (lattice_v->is_flag) {
            std::set<decoding::vertex_t*> affected;
            for (auto dw : gr.get_neighbors(dv)) {
                if (dw->id == BOUNDARY_INDEX)   continue;
                uint dj = from_rlatt[std::make_pair(dw->id, dec_index)];
                cdet_t cdj = std::make_pair(detector_to_base[dj], color_map[dj]);
                stv_t* lattice_w = lattice_structure->get_vertex(CDET_TO_ID(cdj));
                if (!lattice_w->is_flag)    affected.insert(dw);
            }
            for (auto dw : affected) {
                updated_edges.insert(std::make_pair(dv, dw));
                updated_edges.insert(std::make_pair(dw, dv));
                for (auto du : affected) {
                    if (dw == du)   continue;
                    updated_edges.insert(std::make_pair(dw, du));
                }
            }
        } else {
            restricted_syndrome[ddi] = 1;
            vertices_in_syndrome.push_back(dv);
        }
    }

    dec->override_weights.clear();

    for (uint i = 0; i < vertices_in_syndrome.size(); i++) {
        auto dv = vertices_in_syndrome[i];
        for (uint j = i+1; j < vertices_in_syndrome.size(); j++) {
            auto dw = vertices_in_syndrome[j]; 
            auto error_data = gr.get_error_chain_data(dv, dw);
            auto error_chain = error_data.error_chain;

            bool any_modifications = false;
            fp_t updated_weight = 0.0;
            for (uint k = 1; k < error_chain.size(); k++) {
                auto dx = error_chain[k-1];
                auto dy = error_chain[k];
                auto e = gr.get_edge(dx, dy);
                if (updated_edges.count(std::make_pair(dx, dy))) {
                    any_modifications = true;
                    updated_weight += e->edge_weight * 0.5;
                } else {
                    updated_weight += e->edge_weight;
                }
            }
            if (any_modifications) {
                dec->override_weights[std::make_pair(dv, dw)] = updated_weight;
                dec->override_weights[std::make_pair(dw, dv)] = updated_weight;
            }
        }
    }

    auto res = dec->decode_error(restricted_syndrome);
    // Convert the error assignments to the original circuit's detectors.
#ifdef DEBUG
    std::cout << "\tAssignments on lattice restricting color " 
        << c2i(restricted_color) << ":\n";
#endif
    // First, filter out any assignments that go through the boundary.
    for (auto& aa : res.error_assignments) {
        uint& d1 = std::get<0>(aa);
        uint& d2 = std::get<1>(aa);
        if (d1 != BOUNDARY_INDEX)   d1 = from_rlatt[std::make_pair(d1, dec_index)];
        if (d2 != BOUNDARY_INDEX)   d2 = from_rlatt[std::make_pair(d2, dec_index)];
#ifdef DEBUG
        std::cout << "\t\t" << d1 << " <----> " << d2 << "\n";
#endif
    }
    return res;
}

bool
RestrictionDecoder::is_flag(cdet_t cd) {
    cd.first = detector_to_base[cd.first];
    stv_t* v = lattice_structure->get_vertex(CDET_TO_ID(cd));
    return v != nullptr && v->is_flag;
}

std::set<cdetpair_t>
RestrictionDecoder::get_all_edges_in_component(const component_t& comp) {
    std::map<cdet_t, cdet_t> flag_table;

    std::set<cdetpair_t> edges;
    __COLOR cc_color = comp.second;
    for (match_t mate : comp.first) {
        cdet_t cdx = std::get<0>(mate);
        cdet_t cdy = std::get<1>(mate);
#ifdef DEBUG
        std::cout << "\t\t" << cdx.first << "(" << c2i(cdx.second) << ") <----> "
            << cdy.first << "(" << c2i(cdy.second) << ")\n";
#endif
        __COLOR restricted_color = std::get<2>(mate);
        if (restricted_color == cc_color)  {
            continue;
        }

        if (detector_to_base[cdx.first] == detector_to_base[cdy.first]) continue;

        if (is_flag(cdx))   std::swap(cdx, cdy);
        if (is_flag(cdy) && !is_flag(cdx)) {
            if (flag_table.count(cdy)) {
                cdy = flag_table[cdy];
                flag_table.erase(cdy);
#ifdef DEBUG
                std::cout << "\t\t\tnow pairing " << "\t\t" 
                    << cdx.first << "(" << c2i(cdx.second) << ") <----> "
                    << cdy.first << "(" << c2i(cdy.second) << ")\n";
#endif
            } else {
                flag_table[cdy] = cdx;
                continue;
            }
        }

        uint dec_index = 2 - c2i(restricted_color);
        DecodingGraph& gr = rlatt_dec[dec_index]->decoding_graph;
        uint64_t ddx = cdx.first == BOUNDARY_INDEX ? cdx.first : to_rlatt[cdx.first][dec_index];
        uint64_t ddy = cdy.first == BOUNDARY_INDEX ? cdy.first : to_rlatt[cdy.first][dec_index];

        auto vx = gr.get_vertex(ddx);
        auto vy = gr.get_vertex(ddy);
        auto error_chain = gr.get_error_chain_data(vx, vy).error_chain;
        if (error_chain.size() == 0)    continue;
        if (error_chain[0] == vy)   std::reverse(error_chain.begin(), error_chain.end());
        // Now, we must convert the error chain to a path of detectors.
#ifdef DEBUG
        std::cout << "\t\t\tPath:";
#endif
        __COLOR prev_boundary_color = __COLOR::none;
        for (uint dec_index = 1; dec_index < error_chain.size(); dec_index++) {
            auto v1 = error_chain[dec_index-1];
            auto v2 = error_chain[dec_index];

            uint64_t d1 = v1->id;
            uint64_t d2 = v2->id;
            if (d1 == BOUNDARY_INDEX)   std::swap(d1, d2);
            d1 = detector_to_base[from_rlatt[std::make_pair(d1, dec_index)]];
            __COLOR c1 = color_map[d1];
            cdet_t cd1 = std::make_pair(d1, c1);
            __COLOR c2;
            if (d2 == BOUNDARY_INDEX) {
                c2 = get_remaining_color(c1, restricted_color);
                cdet_t cb1 = std::make_pair(d2, c2);
                if (prev_boundary_color != __COLOR::none && prev_boundary_color != c2) {
                    cdet_t tmp = std::make_pair(BOUNDARY_INDEX, prev_boundary_color);
                    insert_pair_into(edges, std::make_pair(cb1, tmp));
                }
                prev_boundary_color = c2;
                // Check if d1 is even adjacent to this boundary.
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
#ifdef DEBUG
                    std::cout << " [ " << cd1.first << "(" << c2i(cd1.second) << ") , "
                            << cb2.first << "(" << c2i(cb2.second) << ") , "
                            << cb1.first << "(" << c2i(cb1.second) << ") ]";
#endif
                    xor_pair_into(edges, std::make_pair(cd1, cb2));
                    xor_pair_into(edges, std::make_pair(cb2, cb1));
                    continue;
                }
            } else {
                d2 = detector_to_base[from_rlatt[std::make_pair(d2, dec_index)]];
                c2 = color_map[d2];
            }
            cdet_t cd2 = std::make_pair(d2, c2);
            cdetpair_t e = std::make_pair(cd1, cd2);
#ifdef DEBUG
            std::cout << " [ " << cd1.first << "(" << c2i(cd1.second) << ") , "
                    << cd2.first << "(" << c2i(cd2.second) << ") ]";
#endif
            xor_pair_into(edges, e);
        }
#ifdef DEBUG
        std::cout << "\n";
#endif
    }
    return edges;
}

std::vector<cdet_t>
RestrictionDecoder::get_detectors_between(cdet_t cd1, cdet_t cd2) {
    auto v1 = lattice_structure->get_vertex(CDET_TO_ID(cd1));
    auto v2 = lattice_structure->get_vertex(CDET_TO_ID(cd2));
    auto path = lattice_matrix[v1][v2];
    std::vector<cdet_t> path_w_cdet;
    for (auto x : path) path_w_cdet.push_back(x->detector);
    return path_w_cdet;
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
    stv_t* v = lattice_structure->get_vertex(CDET_TO_ID(cd));
    for (auto w : lattice_structure->get_neighbors(v)) {
        if (w->is_flag) continue;
        // Filter out neighbors that do not affect the logical state.
        cdet_t cw = w->detector;
        __COLOR restricted_color = get_remaining_color(cd.second, cw.second);
        int dec_index = 2 - c2i(restricted_color);
        uint64_t ddv = cd.first == BOUNDARY_INDEX ? BOUNDARY_INDEX : to_rlatt[cd.first][dec_index];
        uint64_t ddw = cw.first == BOUNDARY_INDEX ? BOUNDARY_INDEX : to_rlatt[cw.first][dec_index];
        DecodingGraph& gr = rlatt_dec[dec_index]->decoding_graph;
        auto vv = gr.get_vertex(ddv);
        auto ww = gr.get_vertex(ddw);
        if (gr.get_error_chain_data(vv, ww).frame_changes.size()
            || ddw == BOUNDARY_INDEX) 
        {
            neighbors.insert(w->detector);
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
    std::cout << "\t\t\t\tApplying corrections on face " 
                << x.first << "(" << c2i(x.second) << "), "
                << y.first << "(" << c2i(y.second) << "), "
                << z.first << "(" << c2i(z.second) << "):";
#endif
    for (uint dec_index = 0; dec_index < detectors.size(); dec_index++) {
        cdet_t cd1 = detectors[dec_index];
        for (uint j = dec_index+1; j < detectors.size(); j++) {
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
                corr[fr] |= 1;
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
