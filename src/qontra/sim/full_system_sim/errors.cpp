/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#include "qontra/sim/full_system_sim.h"

namespace qontra {

void
FullSystemSimulator::inject_timing_error() {
    recalibrate_timing();

    std::vector<uint64_t> qubits;
    std::vector<fp_t> xy_array, z_array;
    for (uint64_t q = 0; q < n_qubits; q++) {
        fp_t t1 = config.timing.t1[q],
            t2 = config.timing.t2[q];
        
        fp_t e_ad = 0.25 * (1 - exp(-elapsed_time/t1));
        fp_t e_pd = 0.5 * (1 - exp(-elapsed_time/t2));

        xy_array.push_back(e_ad);
        z_array.push_back(e_pd-e_ad);

        if (is_recording_stim_instructions) {
            uint32_t _q = static_cast<uint32_t>(q);
            sample_circuit.safe_append_ua("X_ERROR", {_q}, e_ad);
            sample_circuit.safe_append_ua("Y_ERROR", {_q}, e_ad);
            sample_circuit.safe_append_ua("Z_ERROR", {_q}, e_pd-e_ad);
        }
        qubits.push_back(q);
    }
    base_sim->error_channel<&StateSimulator::eX>(qubits, xy_array);
    base_sim->error_channel<&StateSimulator::eY>(qubits, xy_array);
    base_sim->error_channel<&StateSimulator::eZ>(qubits, z_array);
    // If there are any trials with time deltas, handle them now.
    for (auto pair : shot_time_delta_map) {
        uint64_t t = pair.first;
        fp_t delta = pair.second;
        for (uint64_t q = 0; q < n_qubits; q++) {
            fp_t t1 = config.timing.t1[q];
            fp_t t2 = config.timing.t2[q];
            
            fp_t e_ad = 0.25 * (1 - exp(-delta/t1));
            fp_t e_pd = 0.5 * (1 - exp(-delta/t2));

            if (base_sim->get_probability_sample_from_rng() < e_ad)         base_sim->eX(q, t);
            if (base_sim->get_probability_sample_from_rng() < e_ad)         base_sim->eY(q, t);
            if (base_sim->get_probability_sample_from_rng() < e_pd-e_ad)    base_sim->eZ(q, t);
        }
    }
    elapsed_time = 0;
    shot_time_delta_map.clear();
}

void
FullSystemSimulator::inject_idling_error_positive(std::vector<uint64_t> on_qubits, int64_t trial) {
    // Do NOT record the error if trial >= 0.
    std::vector<fp_t> error_rates;
    for (uint64_t q : on_qubits) {
        fp_t e = config.errors.idling[q];
        error_rates.push_back(e);
        if (is_recording_stim_instructions && trial < 0) {
            uint32_t _q = static_cast<uint32_t>(q);
            sample_circuit.safe_append_ua("DEPOLARIZE1", {_q}, e);
        }
    }

    if (trial >= 0) {
        // We are only injecting the error on a single trial.
        for (size_t i = 0; i < error_rates.size(); i++) {
            uint64_t q = on_qubits[i];
            if (base_sim->get_probability_sample_from_rng() < error_rates[i]) base_sim->eDP1(q, trial);
        }
    } else {
        base_sim->error_channel<&StateSimulator::eDP1>(on_qubits, error_rates);
    }
}

void
FullSystemSimulator::inject_idling_error_negative(std::vector<uint64_t> not_on_qubits, int64_t trial) {
    std::vector<uint64_t> on_qubits;
    for (uint64_t i = 0; i < n_qubits; i++) {
        if (std::find(not_on_qubits.begin(), not_on_qubits.end(), i) == not_on_qubits.end()) {
            on_qubits.push_back(i);
        }
    }
    inject_idling_error_positive(on_qubits, trial);
}

}   // qontra
