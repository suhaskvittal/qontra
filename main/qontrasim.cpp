/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#include <qontra/ext/qes.h>
#include <qontra/sim/base/frame_sim.h>
#include <qontra/sim/full_system_sim.h>

#include <vtils/cmd_parse.h>
#include <vtils/filesystem.h>
#include <vtils/ini_parse.h>
#include <vtils/timer.h>

using namespace qontra;
using namespace vtils;

int main(int argc, char* argv[]) {
    std::string help = "usage: ./qontrasim <qes-file>"
                        "--shots <shots> "
                        "--trace-out <trace-folder> "
                        "--stim-out <stim-file> "
                        "--data-out <data-file>\n"
                        "Optional:\n"
                        "\t--subroutines <ini-file>\n";
    CmdParser pp(argc, argv, 1);
    pp.help = help;

    std::string qes_file(argv[1]);

    uint64_t shots;
    std::string trace_folder,
                stim_file,
                data_file,
                ini_file;
    pp.get_uint64("shots", shots, true);
    pp.get_string("trace-out", trace_folder, true);
    pp.get_string("stim-out", stim_file, true);
    pp.get_string("data-out", data_file, true);

    // Initialize simulator:
    FullSystemSimulator sim;

    if (pp.get_string("subroutines", ini_file)) {
        IniParser ini(ini_file);
        const auto& ini_map = ini.get_ini_map();
        for (const auto& p : ini_map.at("__ANON__")) {
            qes::Program<> subroutine = qes::from_file(p.second);
            sim.load_subroutine(p.first, subroutine);
        }
    }

    if (!file_exists(trace_folder)) safe_create_directory(trace_folder);
    if (!file_exists(stim_file)) safe_create_directory(get_parent_directory(stim_file.c_str()));
    if (!file_exists(data_file)) safe_create_directory(get_parent_directory(data_file.c_str()));

    qes::Program<> main_program = qes::from_file(qes_file);
    histogram_t<uint64_t> shots_hist = sim.run_program<FrameSimulator>(main_program, shots);
    histogram_t<double> norm_hist = histogram_normalize(shots_hist);

    std::cout << "Shots Histogram ------------------------\n"
            << shots_hist << "\n"
            << "Probability Histogram ------------------\n"
            << norm_hist << std::endl;
}
