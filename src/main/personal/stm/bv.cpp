/*
 *  author: Suhas Vittal
 *  date:   12 September 2023
 * */

#include "experiments.h"
#include "instruction.h"
#include "parsing/cmd.h"
#include "sim/clifford_sim.h"
#include "sim/manager.h"
#include "tables.h"

using namespace qontra;

int main(int argc, char* argv[]) {
    CmdParser pp(argc, argv);

    MPI_Init(NULL, NULL);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    tables::ErrorAndTiming noise_model_def;

    // Error rates:
    noise_model_def.e_m1w0 = 1e-2;
    noise_model_def.e_m0w1 = 3e-2;
    noise_model_def.e_g1q = 5e-4;
    noise_model_def.e_g2q = 1e-2;
    // Timing:
    noise_model_def.t_g1q = 40;
    noise_model_def.t_g2q = 200;
    noise_model_def.t_ro = 800;
    noise_model_def.t1 = 300e3;
    noise_model_def.t2 = 200e3;

    std::string bv_asm_file;
    uint64_t true_meas_bits;
    uint64_t shots = 100000;

    std::string sim = "clifford";

    if (!pp.get_string("asm", bv_asm_file)) return 1;
    if (!pp.get_uint64("meas", true_meas_bits)) return 1;
    pp.get_uint64("shots", shots);
    pp.get_string("sim", sim);
    schedule_t prog = schedule_from_file(bv_asm_file);
    const uint n = get_number_of_qubits(prog);

    SimManager mgr(prog);
    tables::populate(n, mgr.params.errors, mgr.params.timing, noise_model_def);

    experiments::G_SHOTS_PER_BATCH = shots < 100'000 ? shots : 100'000;

    if (sim == "clifford") {
        CliffordSimulator* csim = new CliffordSimulator(n, experiments::G_SHOTS_PER_BATCH);
        mgr.sim = csim;
    }
    /*else if (sim == "universal") {
        experiments::G_SHOTS_PER_BATCH = 1;
        UniversalSimulator<2>* usim = new UniversalSimulator<2>(n);
        mgr.sim = usim;
    }*/

    mgr.params.always_inject_timing_errors = true;

    auto res = mgr.evaluate_monte_carlo(shots);
    if (world_rank == 0) {
        std::cout << "Entries before: " << res.probability_histogram.size() << "\n";
        // Filter out measurements.
        vlw_t filter_word((true_meas_bits >> 6)+1, 0xffff'ffff'ffff'ffff);
        filter_word[true_meas_bits >> 6] = (1L << (true_meas_bits % 64))-1;
        for (auto it = res.probability_histogram.begin();
                it != res.probability_histogram.end(); ) 
        {
            if (vlw_compare(it->first, filter_word) > 0) {
                it = res.probability_histogram.erase(it);
            } else {
                it++;
            }
        }
        std::cout << "Entries after: " << res.probability_histogram.size() << "\n";
        print_histogram(std::cout, true_meas_bits, res.probability_histogram, 1e-3);
    }
    delete mgr.sim;

    MPI_Finalize();

    return 0;
}
