/*
 * author:  Suhas Vittal
 * date:    27 August 2023
 * */

#include <instruction.h>
#include <parsing/cmd.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace qontra;

#define T_ASM   0
#define T_STIM  1

int get_file_type(std::string filename) {
    const std::vector<std::string> ext{
        ".asm",
        ".stim"
    };
    for (uint i = 0; i < ext.size(); i++) {
        if (filename.find(ext[i]) != std::string::npos) {
            return i;
        }
    }
    return -1;
}

int main(int argc, char* argv[]) {
    CmdParser pp(argc, argv);

    std::string help = "usage: ./converter --input <F> --output <G>\n"
                        "The following conversions are supported:\n"
                        "\t(1) asm to stim\n";
    std::string input_file;
    std::string output_file;

    if (pp.option_set("h")) {
help_exit:
        std::cout << help;
        return 1;
    }

    if (!pp.get_string("input", input_file))    goto help_exit;
    if (!pp.get_string("output", output_file))  goto help_exit;

    int type1 = get_file_type(input_file), type2 = get_file_type(output_file);

    if (type1 == 0) {
        schedule_t sch = schedule_from_file(input_file);
        uint n = get_number_of_qubits(sch);

        // Get error model.
        tables::ErrorAndTiming et;
        ErrorTable errors;
        TimeTable timing;
        tables::populate(n, errors, timing, et);

        stim::Circuit circ = schedule_to_stim(sch, errors, timing);
        std::ofstream out(output_file);
        out << circ.str();
    }
    
    return 0;
}
