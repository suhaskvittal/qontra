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
