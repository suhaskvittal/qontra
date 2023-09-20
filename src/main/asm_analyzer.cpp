/*
 *  author: Suhas Vittal
 *  date:   16 August 2023
 * */

#include "parsing/cmd.h"
#include "instruction.h"

using namespace qontra;

int main(int argc, char* argv[]) {
    CmdParser pp(argc, argv);

    std::string help = 
        "Usage: ./asm_analyzer --asm <file>\n";
    if (pp.option_set("h")) {
        std::cout << help;
        return -1;
    }
    
    std::string file;
    if (!pp.get_string("asm", file)) {
        std::cout << help;
        return -1;
    }

    schedule_t sch = schedule_from_file(file);
    uint64_t pc = 0;
    for (const auto& inst : sch) {
        std::cout << "(pc = " << pc << ") " << inst.str();
        std::cout << "\t@[ ";
        for (auto x : inst.annotations) {
            std::cout << "\"" << x << "\" ";
        }
        std::cout << "]\n";
        pc++;
    }
}
