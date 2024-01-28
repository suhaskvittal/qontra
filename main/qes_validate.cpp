/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#include <qontra/ext/qes.h>
#include <qontra/isa.h>

using namespace qontra;

int main(int argc, char* argv[]) {
    std::string program_file(argv[1]);

    build_isa_from_file(QONTRA_ISA_FILE);

    std::cout << "ISA:\n";
    for (std::string inst : isa_list_instructions()) {
        std::cout << "\t" << inst << "\n";
    }
    std::cout << "\n";

    qes::Program<> program = qes::from_file(program_file);

    uint64_t pc = 0;
    for (const qes::Instruction<>& inst : program) {
        const isa_data_t& dat = isa_get(inst);
        std::cout << "[ pc = " << pc << " ]\n"
                    << "\tparsing \"" << inst << "\" as:\n"
                    << "\t\tinstruction type = " << static_cast<int>(dat.inst_type) << "\n";
        pc++;
    }
    return 0;
}
