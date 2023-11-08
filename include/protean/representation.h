/*
 *  author: Suhas Vittal
 *  date:   20 October 2023
 * */

#ifndef PROTEAN_REPRESENTATION_h
#define PROTEAN_REPRESENTATION_h

#include "defs.h"

#include <iostream>
#include <map>
#include <vector>

namespace qontra {
namespace protean {

struct css_code_data_t {
    std::vector<uint>   data_qubits;
    std::vector<uint>   parity_qubits;
    std::vector<uint>   flag_qubits;
    
    // Each entry of check schedules should be an array of integers such that the i-th
    // entry is the qubit that will have a CNOT with the parity qubit at time step i.
    // 
    // This should be -1 if no such qubit exists.
    std::map<uint, std::vector<int32_t>>    check_schedules;
    uint schedule_depth;
    
    std::vector<uint>   xparity_list;
    std::vector<uint>   zparity_list;

    std::vector<uint>   zflag_list;
    std::vector<uint>   xflag_list;

    std::vector<std::vector<uint>>  x_obs_list;
    std::vector<std::vector<uint>>  z_obs_list;

    // flag_usage[u][v] corresponds to the flag qubit used by data qubit v during check u.
    std::map<uint, std::vector<uint>>   parity_to_flags;
    TwoLevelMap<uint, uint, uint>       flag_usage;

    std::vector<std::vector<uint>> get_observables_and_products(bool get_x_obs) {
        auto base_obs_list = get_x_obs ? x_obs_list : z_obs_list;
        std::vector<std::vector<uint>> all_obs_prods(base_obs_list);

        for (uint i = 0; i < base_obs_list.size(); i++) {
            for (uint j = i+1; j < base_obs_list.size(); j++) {
                std::set<uint> prod(base_obs_list[i].begin(), base_obs_list[i].end());
                for (uint k = 0; k < base_obs_list[j].size(); k++) {
                    uint q = base_obs_list[j][k];
                    if (prod.count(q))  prod.erase(q);
                    else                prod.insert(q);
                }
                if (prod.size() > base_obs_list[i].size() && prod.size() > base_obs_list[j].size()) {
                    continue;
                }
                all_obs_prods.push_back(std::vector<uint>(prod.begin(), prod.end()));
            }
        }
        return all_obs_prods;
    }

    void print_schedule(std::ostream& out) {
        out << "---------------- X Stabilizers -------------------\n";
        for (uint xp : xparity_list) {
            out << xp << "\t|";
            for (int32_t t : check_schedules[xp]) {
                if (t < 0)  out << "\t_";
                else        out << "\t" << t;
            }
            out << "\n";
        }
        out << "---------------- Z Stabilizers -------------------\n";
        for (uint zp : zparity_list) {
            out << zp << "\t|";
            for (int32_t t : check_schedules[zp]) {
                if (t < 0)  out << "\t_";
                else        out << "\t" << t;
            }
            out << "\n";
        }
    }
};

}   // protean
}   // qontra

#endif  // PROTEAN_REPRESENTATION_h
