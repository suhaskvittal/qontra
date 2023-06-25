/*
 *  author: Suhas Vittal
 *  date:   23 June 2023
 * */

#include "sim/control_sim.h"

#include <iostream>
#include <vector>

using namespace qontra;

int main(int argc, char* argv[]) {
    schedule_t  prog;
    const uint n = 17;
    // Build program.
    std::vector<uint>   all_qubits, data_qubits, x_stab, z_stab;
    for (uint i = 0; i < n; i++)    all_qubits.push_back(i);
    for (uint i = 0; i < 9; i++)    data_qubits.push_back(i);
    for (uint i = 9; i < 13; i++)   x_stab.push_back(i);
    for (uint i = 13; i < n; i++)   z_stab.push_back(i);

    prog.push_back({"reset", all_qubits, {}});
    for (uint r = 0; r < 3; r++) {
        prog.push_back({"h", x_stab, {}});
        prog.push_back({"cx", {9, 6, 10, 4, 11, 8, 1, 13, 7, 14, 5, 15}, {}});
        prog.push_back({"cx", {10, 3, 11, 7, 12, 5, 0, 13, 6, 14, 4, 15}, {}});
        prog.push_back({"cx", {10, 0, 11, 4, 12, 2, 3, 14, 1, 15, 7, 16}, {}});
        prog.push_back({"cx", {9, 3, 10, 1, 11, 5, 4, 14, 2, 15, 8, 16}, {}});
        prog.push_back({"h", x_stab, {}});
        prog.push_back({"mnrc", x_stab, {}});
        prog.push_back({"mrc", z_stab, {}});
        if (r > 0) {
            for (uint i = 0; i < 4; i++) {
                prog.push_back({"event", {(i+5), (i+1)}, {}});
            }
        }
        prog.push_back({"reset", x_stab, {}});
        prog.push_back({"reset", z_stab, {}});
        prog.push_back({"nop", all_qubits, {}});
    }
    prog.push_back({"mrc", data_qubits, {}});
    prog.push_back({"obs", {1, 4, 7}, {}});
    prog.push_back({"savem", {}, {}});
    prog.push_back({"done", {}, {}});
    /*
    const uint n = 2;
    prog.push_back({"h", {0}, {}});
    prog.push_back({"cx", {0, 1}, {}});
    prog.push_back({"mrc", {0, 1}, {}});
    prog.push_back({"obs", {1}, {}});
    prog.push_back({"obs", {2}, {}});
    prog.push_back({"savem", {}, {}});
    prog.push_back({"done", {}, {}});
    */
    /*
    const uint n = 5;
    prog.push_back({"reset", {0, 1, 2, 3, 4}, {}});
    prog.push_back({"mrc", {1, 2, 3}, {}});
    prog.push_back({"obs", {1}, {}});
    prog.push_back({"obs", {2}, {}});
    prog.push_back({"obs", {3}, {}});
    prog.push_back({"savem", {}, {}});
    prog.push_back({"done", {}, {}});
    */

    std::cout << to_text(prog) << "\n";

    experiments::G_USE_MPI = false;
    experiments::G_SHOTS_PER_BATCH = 512;

    ControlSimulator sim(n, nullptr, prog);
    tables::ErrorAndTiming params;
    params.e_ro = 0.5;
    tables::populate(n, sim.params.errors, sim.params.timing, params);

    sim.params.verbose = 1;

    sim.run(8192);
    for (auto pair : sim.prob_histograms) {
        for (uint64_t x : pair.first)   std::cout << x << " ";
        std::cout << "\t" << pair.second << "\n";
    }
}

