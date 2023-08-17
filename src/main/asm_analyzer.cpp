/*
 *  author: Suhas Vittal
 *  date:   16 August 2023
 * */

#include "parsing/cmd.h"
#include "instruction.h"

using namespace qontra;

int main(int argc, char* argv[]) {
    CmdParser pp(argc, arv);

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
        std::cout << "(pc = " << pc << ")\t" <<  inst.name;
        for (auto x : inst.operands) {
            std::cout <<  " " << x;
        }
        std::cout << "\n\tannotations:\n";
        for (auto x : inst.annotations) {
            if (x == Annotation::no_errors) {
                std::cout << "\t\tno errors\n";
            } else if (x == Annotation::inject_timing_errors) {
                std::cout << "\t\tinject timing errors\n";
            } else if (x == Annotation::round_start) {
                std::cout << "\t\tround start\n";
            } else if (x == Annotaiton::flag) {
                std::cout << "\t\tflag\n";
            }
        }
    }
}
