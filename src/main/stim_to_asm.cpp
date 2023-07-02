/*
 *  author: Suhas Vittal
 *  date:   30 June 2023
 * */

#include "instruction.h"
#include "parsing/cmd.h"

#include <stim.h>

#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include <stdio.h>

using namespace qontra;

int main(int argc, char* argv[]) {
    CmdParser parser(argc, argv);

    schedule_t prog;

    stim::Circuit circuit;

    std::string help = "Arguments list:\n";
    help += "\tSpecify an output file with --output.\n";
    help += "\tEither pass in a Stim file (--stim) or the code distance (--d)";
    help += " and number of rounds (--r).\n";
    help += "\tPass in (-x) or (-z) to declare the memory experiment.\n";
    help += "\tPass in (-both) to create a memory experiment including events";
    help += " for both X and Z stabilizers.\n";
    help += "\tError rates do not matter.\n";

    bool help_requested = parser.option_set("h");
    if (help_requested) {
help_exit:
        std::cout << help;
        return 0;
    }

    std::string output_file;
    if (!parser.get_string("output", output_file))  goto help_exit;
    
    bool both_stabilizers = parser.option_set("both");
    
    if (parser.option_set("stim")) {
        std::string stim_file;
        parser.get_string("stim", stim_file);
        FILE* fin = fopen(stim_file.c_str(), "r");
        circuit = stim::Circuit::from_file(fin);
        fclose(fin);
    } else {
        // The error rates don't matter.
        uint d, r;
        if (!parser.get_uint32("d", d)) goto help_exit;
        if (!parser.get_uint32("r", r)) goto help_exit;
        bool is_memory_x = parser.option_set("x");
        std::string task = is_memory_x ? "rotated_memory_x" : "rotated_memory_z";
        stim::CircuitGenParameters circ_params(r, d, task);
        circ_params.both_stabilizers = both_stabilizers;
        circuit = generate_surface_code_circuit(circ_params).circuit;
    }
    from_stim_circuit(circuit, prog);
    prog.push_back({"done", {}, {}});
    // Now write the new program to the output file
    std::ofstream out(output_file);
    out << schedule_to_text(prog);
}

