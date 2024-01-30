/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#include <qontra/ext/qes.h>
#include <qontra/sim/base/clifford_sim.h>
#include <qontra/sim/base/frame_sim.h>
#include <qontra/sim/full_system_sim.h>

#include <vtils/cmd_parse.h>
#include <vtils/filesystem.h>
#include <vtils/ini_parse.h>
#include <vtils/timer.h>

using namespace qontra;
using namespace vtils;

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::string help = "usage: ./qontrasim <qes-file> "
                        "--shots <shots> "
                        "--trace-out <trace-folder> "
                        "--stim-out <stim-file> "
                        "--data-out <data-file>\n"
                        "Optional:\n"
                        "\t--config <ini-file>";
    CmdParser pp(argc, argv, 1);
    pp.help = help;
    if (pp.option_set("h")) {
        std::cerr << help << std::endl;
        return 1;
    }
    configure_optimal_batch_size();

    // Initialize simulator:
    FullSystemSimulator sim;

    // Get command line inputs.
    std::string qes_file(argv[1]);

    uint64_t shots;
    std::string ini_file;
    pp.get("shots", shots, true);
    pp.get("trace-out", sim.config.syndrome_output_folder, true);
    pp.get("stim-out", sim.config.stim_output_file, true);
    pp.get("data-out", sim.config.data_output_file, true);

    fp_t p = 0.0;
    pp.get("p", p);
    if (pp.get("config", ini_file)) {
        std::string parent_dir = get_parent_directory(ini_file);

        IniParser ini(ini_file);
        const auto& ini_map = ini.get_ini_map();
        std::cout << "ini keys:";
        for (const auto& p : ini_map) std::cout << " " << p.first;
        std::cout << "\n";
        // First load in the subroutines.
        for (const auto& p : ini_map.at("__ANON__")) {
            qes::Program<> subroutine = qes::from_file(parent_dir + "/" + p.second);
            sim.load_subroutine(p.first, subroutine);
        }
        // Now load in error model data.
        ini.get("Error Model", "p", p);
    }
    qes::Program<> main_program = qes::from_file(qes_file);

    // Setup error model.
    tables::ErrorAndTiming et;
    et *= p*1000;

    const uint64_t n = get_number_of_qubits(main_program);
    tables::populate(n, sim.config.errors, sim.config.timing, et);

    histogram_t<uint64_t> shots_hist = sim.run_program<CliffordSimulator>(main_program, shots);
    histogram_t<double> norm_hist = histogram_normalize(shots_hist);

    if (world_rank == 0) {
        std::cout << "Shots Histogram ------------------------\n"
                << shots_hist << "\n"
                << "Probability Histogram ------------------\n"
                << norm_hist << std::endl;
    }
    MPI_Finalize();
    return 0;
}
