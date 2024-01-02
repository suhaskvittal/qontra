/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include "protean/network.h"

#include <fstream>
#include <iostream>

namespace qontra {
namespace protean {

using namespace net;

void
write_schedule_file(std::string output_file, PhysicalNetwork& network) {
    std::ofstream fout(output_file);
    schedule_t sch = network.make_schedule();
    fout << schedule_to_text(sch) << "\n";
}

void
write_coupling_file(std::string output_file, PhysicalNetwork& network) {
    std::ofstream fout(output_file);
    fout << "Qubit 1,Qubit 2,Processor Layer\n";
    for (sptr<phys_edge_t> pe : network.get_edges()) {
        sptr<phys_vertex_t> pv = network.get_source(pe),
                            pw = network.get_target(pe);
        const size_t layer = pe->tsv_layer;
        fout << pv->id << "," << pw->id << "," << layer << "\n";
    }
}

void
write_role_file(std::string output_file, PhysicalNetwork& network) {
    std::ofstream fout(output_file);
    fout << "Physical Qubit,Roles\n";
    for (sptr<phys_vertex_t> pv : network.get_vertices()) {
        fout << pv->id << ",\"";
        bool first = true;
        for (sptr<raw_vertex_t> rv : pv->role_set) {
            if (!first) fout << " ";
            first = false;

            std::string role;
            if (rv->qubit_type == raw_vertex_t::type::data) {
                role += "d";
            } else if (rv->qubit_type == raw_vertex_t::type::xparity) {
                role += "x";
            } else if (rv->qubit_type == raw_vertex_t::type::zparity) {
                role += "z";
            } else if (rv->qubit_type == raw_vertex_t::type::flag) {
                role += "f";
            } else {
                role += "pr";
            }
            role += std::to_string(rv->id);
            fout << role;
        }
        fout << "\"\n";
    }
}

void
write_tanner_graph_file(std::string output_file, PhysicalNetwork& network) {
    std::ofstream fout(output_file);
    
    using namespace graph;
    RawNetwork raw_net = network.get_raw_connection_network();
    TannerGraph& tanner_graph = raw_net.tanner_graph;
    // First write checks.
    fout << "Parity Role/Logical Observable,Role Operator String,Physical Operator String\n";
    for (sptr<tanner::vertex_t> tpq : tanner_graph.get_checks()) {
        sptr<raw_vertex_t> rpq = raw_net.v_tanner_raw_map.at(tpq);
        if (rpq->qubit_type == raw_vertex_t::type::xparity) {
            fout << "x";
        } else {
            fout << "z";
        }
        fout << rpq->id;
        std::string role_check_string, phys_check_string;
        bool first = true;
        for (sptr<tanner::vertex_t> tdq : tanner_graph.get_neighbors(tpq)) {
            sptr<raw_vertex_t> rdq = raw_net.v_tanner_raw_map.at(tdq);
            sptr<phys_vertex_t> pdq = network.role_to_phys[rdq];
            
            if (!first) {
                role_check_string += ",";
                phys_check_string += ",";
            }
            first = false;
            role_check_string += "d" + std::to_string(rdq->id);
            phys_check_string += std::to_string(pdq->id);
        }
        fout << ",\"" << role_check_string << "\",\"" << phys_check_string << "\"\n";
    }
    // Now write logical observables.
    for (size_t c = 0; c <= 1; c++) {
        bool is_x_obs = (c == 0);
        size_t k = 0;
        for (auto obs : tanner_graph.get_obs(is_x_obs)) {
            fout << "o" << (is_x_obs ? "x" : "z") << k;

            std::string role_obs_string, phys_obs_string;
            bool first = true;
            for (sptr<tanner::vertex_t> tv : obs) {
                if (!first) {
                    role_obs_string += ",";
                    phys_obs_string += ",";
                }
                first = false;
                sptr<raw_vertex_t> rv = raw_net.v_tanner_raw_map.at(tv);
                sptr<phys_vertex_t> pv = network.role_to_phys[rv];
                role_obs_string += "d" + std::to_string(rv->id);
                phys_obs_string += std::to_string(pv->id);
            }
            fout << ",\"" << role_obs_string << "\",\"" << phys_obs_string << "\"\n";
            k++;
        }
    }
}

void
write_flag_assignment_file(std::string output_file, PhysicalNetwork& network) {
    std::ofstream fout(output_file);
    fout << "Parity Role,D:D:F Role List\n";

    RawNetwork raw_net = network.get_raw_connection_network();
    for (auto& p1 : raw_net.flag_assignment_map) {
        sptr<raw_vertex_t> rpq = p1.first;
        if (rpq->qubit_type == raw_vertex_t::type::xparity) {
            fout << "x";
        } else {
            fout << "z";
        }
        fout << rpq->id << ",";

        std::map<sptr<raw_vertex_t>, sptr<raw_vertex_t>>
            flag_to_first_data;
        std::string flag_string;
        bool first = true;
        for (auto& p2 : p1.second) {
            sptr<raw_vertex_t> rdq1 = p2.first,
                               rfq = p2.second;
            if (!flag_to_first_data.count(rfq)) {
                flag_to_first_data[rfq] = rdq1;
            } else {
                // The pair should be complete: write it.
                if (!first) flag_string += ",";
                first = false;

                sptr<raw_vertex_t> rdq2 = flag_to_first_data[rfq];
                flag_string += "d" + std::to_string(rdq1->id) + ":"
                                + "d" + std::to_string(rdq2->id) + ":"
                                + "f" + std::to_string(rfq->id);
            }
        }
        fout << "\"" + flag_string + "\"\n";
    }
}


}   // protean
}   // qontra
