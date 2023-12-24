
#include <instruction.h>
#include <protean/scheduler.h>
#include <protean/utils.h>

#include <fstream>
#include <iostream>

using namespace qontra;
using namespace protean;

int main() {
    css_code_data_t dd;

    for (uint i = 0; i < 24; i++) dd.data_qubits.push_back(i);
        
    std::vector<int32_t> checks[] = {
        {-1, -1, -1, -1, -1, -1, -1, -1, 3, 4, 2, 0, 1, 5},
        {-1, -1, -1, -1, -1, -1, -1, -1, 1, 2, 0, 5, 3, 4},
        {-1, -1, -1, -1, -1, -1, -1, -1, 20, 19, 21, 23, 22, 18},
        {-1, -1, -1, -1, -1, -1, -1, -1, 22, 21, 23, 18, 20, 19},
        {-1, -1, -1, -1, -1, -1, -1, -1, 8, 11, 9, 7, 6, 10},
        {-1, -1, -1, -1, -1, -1, -1, -1, 6, 9, 7, 10, 8, 11},
        {-1, -1, -1, -1, -1, -1, -1, -1, 15, 12, 14, 16, 17, 13},
        {-1, -1, -1, -1, -1, -1, -1, -1, 17, 14, 16, 13, 15, 12},

        {22, 3, 7, 9, 20, 16, 1, 14, -1, -1, -1, -1, -1, -1},
        {3, 14, 1, 16, 9, 22, 20, 7, -1, -1, -1, -1, -1, -1},
        {17, 8, 23, 21, 15, 0, 6, 2, -1, -1, -1, -1, -1, -1},
        {8, 2, 6, 0, 21, 17, 15, 23, -1, -1, -1, -1, -1, -1},
        {9, 20, 12, 5, 3, 11, 14, 18, -1, -1, -1, -1, -1, -1},
        {20, 18, 14, 11, 5, 9, 3, 12, -1, -1, -1, -1, -1, -1},
        {21, 15, 4, 10, 8, 19, 2, 13, -1, -1, -1, -1, -1, -1},
        {15, 13, 2, 19, 10, 21, 8, 4, -1, -1, -1, -1, -1, -1},
    };

    std::vector<uint> obs_list[] = {
        {20, 22, 7, 9},
        {15, 17, 21, 23},
        {0, 2, 6, 8},
        {1, 3, 14, 16},
        {2, 4, 8, 10},
        {3, 5, 12, 14},
        {9, 11, 18, 20},
        {13, 15, 19, 21}
    };

    uint k = 0;
    for (uint i = 24; i < 40; i++) {
        if (i % 2 == 0) dd.zparity_list.push_back(i);
        else            dd.xparity_list.push_back(i);
        dd.parity_qubits.push_back(i);
        dd.check_schedules[i] = checks[k++];
    }
    dd.schedule_depth = 14;

    for (auto obs : obs_list) {
        dd.x_obs_list.push_back(obs);
        dd.z_obs_list.push_back(obs);
    }

    dd.print_schedule(std::cout);
    dd = protean::make_fault_tolerant_simple(dd);

    schedule_t mxp = write_memory_experiment(dd, 4, false);

    std::ofstream fout("test.csv");
    fout << schedule_to_text(mxp);
}
