/*
 *  author: Suhas Vittal
 *  date:   23 January 2024
 * */

namespace qontra {

inline void
FullSystemSimulator::execute_routine(const qes::Program<>& program) {
    program_status_t status;
    while (status.pc < program.size()) {
        const qes::Instruction<>& instruction = program.at(status.pc);
        // Check if the instruction is requesting the execution of the
        // microcode.
        if (instruction.get_name() == "call") {
            // Switch into the subroutine.
            std::string subroutine_name = instruction.get<std::string>(1);
            execute_routine(subroutine_map.at(subroutine_name));
            status.pc++;
        } else {
            read_next_instruction(current_program, main_status);
            if (status.return_if_waiting_trials.popcnt() == current_shots) break;
        }
    }
}

}   // qontra
