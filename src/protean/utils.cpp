/*
 *  author: Suhas Vittal
 *  date:   26 October 2023
 * */

#include "protean/utils.h"

namespace qontra {
namespace protean {

schedule_t
write_memory_experiment(css_code_data_t code_data, uint rounds, bool is_memory_x) {
    schedule_t prog;
    //
    // PROLOGUE.
    //
    std::vector<uint> all_qubits(code_data.data_qubits);
    for (auto q : code_data.parity_qubits) all_qubits.push_back(q);
    for (auto q : code_data.flag_qubits) all_qubits.push_back(q);

    prog.push_back(Instruction::gate("reset", all_qubits));
    if (is_memory_x) {
        prog.push_back(Instruction::gate("h", code_data.data_qubits));
    }
    //
    // BODY
    //
    std::vector<uint> all_flag_parity_qubits;
    std::map<uint, uint> qubit_to_meas_time;
    uint64_t meas_per_round = 0;
    for (uint q : code_data.zparity_list) {
        all_flag_parity_qubits.push_back(q);
        qubit_to_meas_time[q] = meas_per_round++;
    }
    for (uint q : code_data.xparity_list) {
        all_flag_parity_qubits.push_back(q);
        qubit_to_meas_time[q] = meas_per_round++;
    }
    for (uint q : code_data.zflag_list) {
        all_flag_parity_qubits.push_back(q);
        qubit_to_meas_time[q] = meas_per_round++;
    }
    for (uint q : code_data.xflag_list) {
        all_flag_parity_qubits.push_back(q);
        qubit_to_meas_time[q] = meas_per_round++;
    }
    Instruction flagparitymeas = Instruction::gate("measure", all_flag_parity_qubits);
    Instruction flagparityreset = Instruction::gate("reset", all_flag_parity_qubits);
    // Now add all the instructions to the schedule.
    std::vector<uint> event_qubits;
    if (is_memory_x) {
        for (uint q : code_data.xparity_list) event_qubits.push_back(q);
        for (uint q : code_data.zflag_list) event_qubits.push_back(q);
    } else {
        for (uint q : code_data.zparity_list) event_qubits.push_back(q);
        for (uint q : code_data.xflag_list) event_qubits.push_back(q);
    }
    uint64_t event_ctr = 0;
    schedule_t body = write_syndrome_extraction_ops(code_data);
    for (uint r = 0; r < rounds; r++) {
        for (auto inst : body) prog.push_back(inst);
        prog.push_back(flagparitymeas);
        prog.push_back(flagparityreset);
        // Create detection events.
        for (uint i = 0; i < event_qubits.size(); i++) {
            uint q = event_qubits[i];
            uint mt = qubit_to_meas_time[q];
            if (r > 0) {
                uint mt1 = mt + r*meas_per_round;
                uint mt2 = mt + (r-1)*meas_per_round;
                prog.push_back(Instruction::event(event_ctr++, {mt1, mt2}));
            } else {
                prog.push_back(Instruction::event(event_ctr++, {mt}));
            }
        }
    }
    const uint64_t total_measurements = rounds*meas_per_round;
    //
    // Epilogue
    //
    if (is_memory_x) {
        prog.push_back(Instruction::gate("h", code_data.data_qubits));
    }
    // We need to record when we measure the data qubits when creating the final
    // detection events and the observables.
    for (uint i = 0; i < code_data.data_qubits.size(); i++) {
        uint q = code_data.data_qubits[i];
        qubit_to_meas_time[q] = total_measurements+i;
    }
    prog.push_back(Instruction::gate("measure", code_data.data_qubits));
    // Create final detection events.
    std::vector<uint> relevant_checks;
    if (is_memory_x) relevant_checks = code_data.xparity_list;
    else             relevant_checks = code_data.zparity_list;
    for (uint pq : relevant_checks) {
        std::vector<int32_t> supp = code_data.check_schedules[pq];
        std::vector<uint> mtimes;
        for (int32_t dq : supp) {
            if (dq < 0) continue;
            uint mt = qubit_to_meas_time[dq];
        }
        mtimes.push_back((rounds-1)*meas_per_round + qubit_to_meas_time[pq]);
        prog.push_back(Instruction::event(event_ctr++, mtimes));
    }
    // Create logical observables.
    std::vector<std::vector<uint>> obs_list;
    if (is_memory_x)    obs_list = code_data.x_obs_list;
    else                obs_list = code_data.z_obs_list;
    uint64_t obs_ctr = 0;
    for (auto obs : obs_list) {
        std::vector<uint> mtimes;
        for (uint dq : obs) mtimes.push_back(qubit_to_meas_time[dq]);
        prog.push_back(Instruction::obs(obs_ctr++, mtimes));
    }
    // The schedule is done -- we can return now.
    return prog;
}

schedule_t
write_syndrome_extraction_ops(css_code_data_t code_data) {
    std::vector<uint> hlist(code_data.xparity_list);
    for (auto q : code_data.zflag_list) hlist.push_back(q);

    // These are the H gates used at the start and end of the round.
    Instruction hstart = Instruction::gate("h", hlist);
    hstart.annotations.insert("inject_timing_error");
    Instruction hend = Instruction::gate("h", hlist);
    // This is the set of CX gates required to entangle/disentangle the flags
    // from the parity qubits.
    std::vector<Instruction> flagcx_array;
    int index = 0;
    while (true) {
        std::vector<uint> qarr;
        bool had_flag = false;
        for (auto pq : code_data.zparity_list) {
            std::vector<uint> flags = code_data.parity_to_flags[pq];
            if (index >= flags.size())  continue;
            had_flag = true;

            uint fq = flags[index];
            qarr.push_back(fq);
            qarr.push_back(pq);
        }
        for (auto pq : code_data.xparity_list) {
            std::vector<uint> flags = code_data.parity_to_flags[pq];
            if (index >= flags.size())  continue;
            had_flag = true;

            uint fq = flags[index];
            qarr.push_back(pq);
            qarr.push_back(fq);
        }
        if (!had_flag) break;
        flagcx_array.push_back(Instruction::gate("cx", qarr));
        index++;
    }

    // Finally, the proper syndrome extraction CNOTs.
    std::vector<Instruction> extcx_array;
    index = 0;
    while (true) {
        std::vector<uint> qarr;
        bool had_op = false;
        for (auto pq : code_data.zparity_list) {
            std::vector<int32_t> supp = code_data.check_schedules[pq];
            if (index >= supp.size())   continue;
            had_op = true;
            if (supp[index] < 0) continue;

            uint dq = supp[index];
            // This may have a flag qubit -- check if it does.
            if (code_data.flag_usage[pq].count(dq)) {
                uint fq = code_data.flag_usage[pq][dq];
                qarr.push_back(dq);
                qarr.push_back(fq);
            } else {
                qarr.push_back(dq);
                qarr.push_back(pq);
            }
        }

        for (auto pq : code_data.xparity_list) {
            std::vector<int32_t> supp = code_data.check_schedules[pq];
            if (index >= supp.size())   continue;
            had_op = true;
            if (supp[index] < 0) continue;

            uint dq = supp[index];
            // This may have a flag qubit -- check if it does.
            if (code_data.flag_usage[pq].count(dq)) {
                uint fq = code_data.flag_usage[pq][dq];
                qarr.push_back(fq);
                qarr.push_back(dq);
            } else {
                qarr.push_back(pq);
                qarr.push_back(dq);
            }
        }
        if (!had_op) break;
        extcx_array.push_back(Instruction::gate("cx", qarr));
        index++;
    }
    // Build schedule.
    schedule_t body;
    body.push_back(hstart);
    for (auto cx_inst : flagcx_array)   body.push_back(cx_inst);
    for (auto cx_inst : extcx_array)    body.push_back(cx_inst);
    for (auto cx_inst : flagcx_array)   body.push_back(cx_inst);
    body.push_back(hend);
    return body;
}

}
}    // qontra
