/*
 *  author: Suhas Vittal
 *  date:   16 August 2023
 * */

#include "parsing/cmd.h"
#include "instruction.h"
#include "sim/enumerator.h"

using namespace qontra;

void f_syntax_analysis(const schedule_t& prog) {
    uint64_t pc = 0;
    for (const auto& inst : prog) {
        std::cout << "(pc = " << pc << ") " << inst.str();
        std::cout << "\tannotations = [ ";
        for (auto x : inst.annotations) {
            std::cout << "\"" << x << "\" ";
        }
        std::cout << "]\n";
        std::cout << "\tproperties = [ ";
        for (auto pair : inst.properties) {
            std::cout << "\"" << pair.first << "\"(i=" << pair.second.ival 
                << ", f=" << pair.second.fval << ") ";
        }
        std::cout << "]\n";
        pc++;
    }
}

void f_error_enumeration(const schedule_t& prog) {
    MPI_Init(NULL, NULL);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    auto record_list = enumerate_errors(prog);
    if (world_rank == 0) {
        write_recorded_errors_to(std::cout, record_list);
    }
}

int main(int argc, char* argv[]) {
    CmdParser pp(argc, argv);

    std::string help = 
        "Usage: ./asm_analyzer -mode --asm <file>\n"
        "\tLegal modes: syntax (-s), error (-e)\n";

    if (pp.option_set("h")) {
        std::cout << help;
        return -1;
    }
    
    std::string file;
    if (!pp.get_string("asm", file)) {
        std::cout << help;
        return -1;
    }
    bool syntax_analysis = pp.option_set("s");
    bool error_enumeration = pp.option_set("e");

    schedule_t prog = schedule_from_file(file);

    if (syntax_analysis)        f_syntax_analysis(prog);
    else if (error_enumeration) f_error_enumeration(prog);
}
