/* author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include "qontra/protean/io.h"

#include <fstream>
#include <iostream>

namespace qontra {

using namespace graph;

namespace protean {

using namespace net;

bool
update_network(std::string pass_string, PhysicalNetwork& network, bool verbose) {
    bool network_was_changed = false;
    if (verbose) {
        std::cout << "[ update_network ] running on pass_string = " << pass_string << std::endl;
    }

    size_t nesting_depth = 0;
    bool capture_string_is_pass_group = false;
    bool capture_string_changed_network = false;
    std::string capture_string;

    // This is some jank syntax, but it typedefs pass_t as a function of
    // PhysicalNetwork that takes in void and returns bool.
    typedef bool (PhysicalNetwork::*pass_t)(void);

#define PAIR(x, y) std::make_pair<std::string, pass_t>(x, &y)
    
    const std::map<std::string, pass_t> PASSES{
        PAIR("jid", PhysicalNetwork::join_qubits_with_identical_support),
        PAIR("jpa", PhysicalNetwork::join_qubits_with_partial_support),
        PAIR("fla", PhysicalNetwork::make_flags),
        PAIR("prx", PhysicalNetwork::add_connectivity_reducing_proxies),
        PAIR("con", PhysicalNetwork::contract_small_degree_qubits),
        PAIR("ral", PhysicalNetwork::reallocate_edges),
        PAIR("rcr", PhysicalNetwork::recompute_cycle_role_maps),
        PAIR("rlb", PhysicalNetwork::relabel_qubits)
    };

    auto call_pass = [&] (std::string p)
    {
        if (verbose) {
            std::cout << "[ update_network ] calling pass " << p << std::endl;
        }
        pass_t pass = PASSES.at(p);
        return (network.*pass)();
    };

    for (size_t i = 0; i <= pass_string.size(); i++) {
        char c = i == pass_string.size() ? '.' : pass_string[i];
        // Execute contents of capture string.
        if (nesting_depth == 0 && capture_string.size() >= 3) {
            capture_string_changed_network = false;
            do {
                if (capture_string_is_pass_group) {
                    capture_string_changed_network = update_network(capture_string, network, verbose);
                } else if (PASSES.count(capture_string)) {
                    capture_string_changed_network = call_pass(capture_string);
                }
                network_was_changed |= capture_string_changed_network;
            } while (c == '+' && capture_string_changed_network);
            capture_string.clear();
            capture_string_is_pass_group = false;
        }
        // Now, update capture_string.
        if (c == ')') { nesting_depth--; capture_string_is_pass_group = true; }
        if (nesting_depth > 0) {
            // All characters are appended without regard.
            capture_string.push_back(c);
        } else {
            if (c == '.' || c == ',' || c == ';' || c == ':' || c == '|'
                || c == '(' || c == ')' || c == '+') 
            {
               // Do nothing.
            } else if (c >= 'A' && c <= 'Z') {
                c += 'a' - 'A';
                capture_string.push_back(c);
            } else if (c >= 'a' && c <= 'z') {
                capture_string.push_back(c);
            } else {
                std::cerr << "[ update_network ] unrecognized character \'" << c << "\'" << std::endl;
            }
        }
        if (c == '(') nesting_depth++;
    }
    return network_was_changed;
}

void
write_stats_file(std::string output_file, PhysicalNetwork& network) {
    std::ofstream fout(output_file);
    fout << "Property,Value\n"
            << "Qubits," << network.n() << "\n"
            << "Couplings," << network.m() << "\n"
            << "Mean Degree," << network.get_mean_connectivity() << "\n"
            << "Max Degree," << network.get_max_connectivity() << "\n"
            << "Thickness," << network.get_thickness() << "\n";
    fout.flush();
}

void
write_schedule_file(std::string output_file, PhysicalNetwork& network) {
    qes::Program<> schedule = network.make_schedule();
    qes::to_file(output_file, schedule);
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
    fout << "Physical Qubit,Degree,Roles\n";
    for (sptr<phys_vertex_t> pv : network.get_vertices()) {
        fout << pv->id << "," << network.get_degree(pv) << ",\"";
        bool first = true;
        for (sptr<raw_vertex_t> rv : pv->role_set) {
            if (!first) fout << ",";
            first = false;
            fout << print_v(rv) << "(" << pv->cycle_role_map.at(rv) << ")";
        }
        fout << "\"\n";
    }
    fout.flush();
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
    fout.flush();
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
    fout.flush();
}


}   // protean
}   // qontra
