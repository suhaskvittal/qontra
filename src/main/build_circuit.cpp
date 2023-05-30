/*
 *  author: Suhas Vittal
 *  date:   30 May 2023
 * */

#include "defs.h"
#include "parsing/cmd.h"

#include <stim.h>

#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    qontra::CmdParser parser(argc, argv);

    std::string help = "Arguments list:\n";
    help += "\tMemory experiment type (X (-x) or Z (-z) default is -z)\n";
    help += "\tRounds (--rounds, default is distance)\n";
    help += "\tDistance (--distance)\n";
    help += "\tOutput file (--output-file, should be a .stim file)\n";
    help += "\tBefore round data depolarization mean (--brddm)\n";
    help += "\tBefore round data depolarization std (--brdds, default 0)\n";
    help += "\tAfter 1-qubit clifford depolarization mean (--a1cdm)\n";
    help += "\tAfter 1-qubit clifford deoplarization std (--a1cds, default 0)\n";
    help += "\tAfter 2-qubit clifford depolarization mean (--a2cdm)\n";
    help += "\tAfter 2-qubit clifford depolarization std (--a2cds, default 0)\n";
    help += "\tBefore measurement flip probability mean (--bmfpm)\n";
    help += "\tBefore measurement flip probability std (--bmfps, default 0)\n";
    help += "\tAfter reset flip probability mean (--arfpm)\n";
    help += "\tAfter reset flip probability std (--afrps, default 0)\n";

    bool help_requested = parser.option_set("h");
    if (help_requested) {
help_exit:
        std::cout << help;
        return 0;
    }

    std::string output_file;
    uint distance;
    uint rounds;
    fp_t brddm, brdds, a1cdm, a1cds, a2cdm, a2cds, bmfpm, bmfps, arfpm, arfps;
    if (!parser.get_float("brddm", brddm)
        || !parser.get_float("a1cdm", a1cdm)
        || !parser.get_float("a2cdm", a2cdm)
        || !parser.get_float("bmfpm", bmfpm)
        || !parser.get_float("arfpm", arfpm)
        || !parser.get_uint32("distance", distance)
        || !parser.get_string("output-file", output_file))
    {
        goto help_exit;
    }

    if (!parser.get_float("brdds", brdds))  brdds = 0;
    if (!parser.get_float("a1cds", a1cds))  a1cds = 0;
    if (!parser.get_float("a2cds", a2cds))  a2cds = 0;
    if (!parser.get_float("bmfps", bmfps))  bmfps = 0;
    if (!parser.get_float("arfps", arfps))  arfps = 0;

    if (!parser.get_uint32("rounds", rounds))   rounds = distance;

    bool is_memory_x = parser.option_set("x");
    std::string task = is_memory_x ? "rotated_memory_x" : "rotated_memory_z";
    stim::CircuitGenParameters params(rounds, distance, task);

    params.before_round_data_depolarization = brddm;
    params.after_clifford_sq_depolarization = a2cdm;
    params.after_clifford_depolarization = a2cdm;
    params.before_measure_flip_probability = bmfpm;
    params.after_reset_flip_probability = arfpm;

    params.before_round_data_depolarization_stddev = brdds;
    params.after_clifford_sq_depolarization_stddev = a1cds;
    params.after_clifford_depolarization_stddev = a2cds;
    params.before_measure_flip_probability_stddev = bmfps;
    params.after_reset_flip_probability_stddev = arfps;

    auto gen_circ = stim::generate_surface_code_circuit(params);
    // Print layout to stdout
    std::cout << gen_circ.layout_str() << "\n";
    // Write to stim file.
    std::filesystem::path output_path(output_file);
    std::filesystem::path output_folder = output_path.parent_path();
    qontra::safe_create_directory(output_folder);

    std::ofstream out(output_path);
    out << gen_circ.circuit << "\n";
}
