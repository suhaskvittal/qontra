/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

#include "graph/lattice_prebuilts.h"

graph::LatticeGraph
surface_code_lattice_graph(uint d) {
    LatticeGraph gr;

    // Make qubits.
    const uint n_data = d*d;
    const uint n_parity = d*d-1;
    for (uint i = 0; i < n_data; i++) {
        sptr<vertex_t> v = std::make_shared<vertex_t>();
        v->id = i;
        v->qubit_type = vertex_t::type::data;
        gr.add_vertex(v);
    }

    typedef std::tuple<int, int>  coord_t;
    std::map<coord_t, sptr<vertex_t>> loc_to_check_map;
    // Make boundary checks.
    uint i = n_data;
    for (uint r = 1; r < d; r += 2) {
        auto pv1 = std::make_shared<vertex_t>();
        pv1->id = (i++);
        pv1->qubit_type = vertex_t::type::zparity;
        auto pv2 = std::make_shared<vertex_t>();
        pv2->id = (i++);
        pv2->qubit_type = vertex_t::type::zparity;

        gr.add_vertex(pv1);
        gr.add_vertex(pv2);

        loc_to_check_map[std::make_tuple(r, 0)] = pv1;
        loc_to_check_map[std::make_tuple(r+1, d)] = pv2;
    }

    for (uint c = 1; c < d; c += 2) {
        auto pv1 = std::make_shared<vertex_t>();
        pv1->id = (i++);
        pv1->qubit_type = vertex_t::type::xparity;
        auto pv2 = std::make_shared<vertex_t>();
        pv2->id = (i++);
        pv2->qubit_type = vertex_t::type::xparity;

        gr.add_vertex(pv1);
        gr.add_vertex(pv2);

        loc_to_check_map[std::make_tuple(d, c)] = pv1;
        loc_to_check_map[std::make_tuple(0, c+1)] = pv2;
    }
    // Now do checks in the bulk.
    for (uint r = 1; r < d; r++) {
        for (uint c = 1; c < d; c++) {
            auto pv = std::make_shared<vertex_t>();
            pv->id = (i++);
            pv->qubit_type = (r+c) & 0x1 ? vertex_t::type::zparity : vertex_t::type::xparity;

            gr.add_vertex(pv);
            
            loc_to_check_map[std::make_tuple(r, c)] = pv;
        }
    }
    // Now create edges between data and parity qubits.
    const int64_t _d = d;
    for (auto& pair : loc_to_check_map) {
        coord_t crd = pair.first;
        int r = std::get<0>(crd),
            c = std::get<1>(crd);
        auto pv = pair.second;

        //  q1      q2
        //      C
        //  q3      q4

        // 1, 2 --> 3 4 6 7

        int64_t q1 = _d*(c-1) + r-1,
                q2 = _d*c + r-1,
                q3 = _d*(c-1) + r,
                q4 = _d*c + r;
        if (r == 0) { q1 = -1; q2 = -1; }
        if (r == d) { q3 = -1; q4 = -1; }
        if (c == 0) { q1 = -1; q3 = -1; }
        if (c == d) { q2 = -1; q4 = -1; }

        int64_t order[] = { q4, q2, q3, q1 };
        if (pv->qubit_type == vertex_t::type::xparity) std::swap(order[1], order[2]);
        for (uint j = 0; j < 4; j++) {
            if (order[j] < 0 || order[j] >= n_data) {
                continue;
            }
            auto dv = gr.get_vertex((uint64_t)order[j]);
            auto e = std::make_shared<edge_t>();
            e->src = std::static_pointer_cast<void>(pv);
            e->dst = std::static_pointer_cast<void>(dv);
            e->cx_time = j;
            gr.add_edge(e);
        }
    }
    // Create observables.
    std::vector<sptr<vertex_t>> x_obs, z_obs;
    for (uint j = 0; j < d; j++) {
        x_obs.push_back(gr.get_vertex(j));
        z_obs.push_back(gr.get_vertex(d*j + d-1));
    }
    gr.x_obs_list.push_back(x_obs);
    gr.z_obs_list.push_back(z_obs);
    return gr;
}

