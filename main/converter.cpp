/*
 * author:  Suhas Vittal
 * date:    27 August 2023
 * */

#include <qontra/ext/stim.h>
#include <qontra/ext/qes.h>
#include <vtils/cmd_parse.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace qontra;

int get_file_type(std::string filename) {
    const std::vector<std::string> ext{
        ".qes",
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
    vtils::CmdParser pp(argc, argv);

    std::string help = "usage: ./converter --in <F> --out <G>\n"
                        "The following conversions are supported:\n"
                        "\t(1) qes to stim\n"
                        "\n"
                        "Optional arguments:\n"
                        "\t--p (physical error rate, default 0.001)";
    std::string input_file;
    std::string output_file;

    fp_t p = 1e-3;

    if (pp.option_set("h")) {
help_exit:
        std::cerr << help << std::endl;
        return 1;
    }

    if (!pp.get_string("in", input_file))    goto help_exit;
    if (!pp.get_string("out", output_file))  goto help_exit;

    pp.get_float("p", p);

    int type1 = get_file_type(input_file), type2 = get_file_type(output_file);

    if (type1 == 0) {
        qes::Program<> program = qes::from_file(input_file);
        size_t n = get_number_of_qubits(program);

        // Get error model.
        tables::ErrorAndTiming et;
        ErrorTable errors;
        TimeTable timing;
        tables::populate(n, errors, timing, et);
        et *= 1000*p;

        DetailedStimCircuit circuit = DetailedStimCircuit::from_qes(program, errors, timing);
        std::ofstream out(output_file);
        out << circuit.str();
    }
    
    return 0;
}
