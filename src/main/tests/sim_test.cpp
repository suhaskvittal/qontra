/*
 *  author: Suhas Vittal
 *  date:   19 August 2023
 * */

#include "parsing/cmd.h"
#include "sim/clifford_sim.h"
#include "sim/manager.h"

using namespace qontra;

int main(int argc, char* argv[]) {
    CmdParser pp(argc, argv);

    std::string asm_file;
    pp.get_string("asm", asm_file);

    schedule_t program = schedule_from_file(asm_file);
    uint n_qubits = get_number_of_qubits(program);
    SimManager mgr(program);
    
    CliffordSimulator state_sim(n_qubits, experiments::G_SHOTS_PER_BATCH);
    mgr.sim = &state_sim;

    tables::ErrorAndTiming et;
    tables::populate(n_qubits, mgr.params.errors, mgr.params.timing, et);

    experiments::G_USE_MPI = false;

    auto res = mgr.evaluate_monte_carlo(100'000);
    for (auto pair : res.probability_histogram) {
        vlw_t w = pair.first;
        uint64_t cnt = pair.second;

        for (uint i = 0; i < w.size(); i++) {
            if (i > 0) std::cout << "|";
            std::cout << w[i];
        }
        std::cout << " : " << cnt << "\n";
    }
}
