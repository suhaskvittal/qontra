/*
 *  author: Suhas Vittal
 *  date:   8 August 2022
 * */

#include "lilliput.h"

LILLIPUT::LILLIPUT(const stim::Circuit& circuit, Decoder * decoder_p,
        const LILLIPUTParams& params, std::mt19937_64& rng)
:Decoder(circuit), 
    lut_misses(0),
    params(params), 
    lut(), 
    is_lut_entry_valid(),
    backing_decoder_p(decoder_p),
    lut_memory_overhead(0.0)
{
    backing_decoder_p->match_detectors_less_than = params.detectors_per_round;
    populate_lut(rng);
}

std::string
LILLIPUT::name() {
    return "LILLIPUT";
}

bool 
LILLIPUT::is_software() {
    return false;
}

DecoderShotResult
LILLIPUT::decode_error(const std::vector<uint8_t>& syndrome) {
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();
    // Need to perform sliding window and state updates.
    uint window_size = backing_decoder_p->circuit.count_detectors();

    std::vector<uint8_t> correction(n_observables, 0);
    std::vector<uint8_t> state(params.detectors_per_round, 0);

    // Update syndrome as we neutralize errors.
    std::vector<uint8_t> running_syndrome(syndrome);

    uint start = 0;
    uint end = window_size;

    std::map<uint, uint> matching;
    while (end < syndrome.size()) {
        std::vector<uint8_t> sub_syndrome(
                running_syndrome.begin() + start, running_syndrome.begin() + end);
        // Apply state
        for (uint i = 0; i < params.detectors_per_round; i++) {
            sub_syndrome[i] ^= state[i];
        }
        // First check if oldest round has any flipped bits.
        bool skip_round = true;
        for (uint i = 0; i < params.detectors_per_round; i++) {
            skip_round = skip_round && !sub_syndrome[i];
        }
        if (skip_round) {
            start += params.detectors_per_round;
            end += params.detectors_per_round;
            continue;
        }
        // Otherwise, decode this window and update state.
        auto res = get_lut_entry(sub_syndrome);
        if (res.correction.size() == 0) {
            start += params.detectors_per_round;
            end += params.detectors_per_round;
            continue;
        }
        std::set<uint> visited;
        for (auto pair : res.matching) {
            if (visited.count(pair.first) || visited.count(pair.second)) {
                continue;
            }
            if (pair.first >= params.detectors_per_round) continue;
                
            if (pair.first != BOUNDARY_INDEX) {
                running_syndrome[pair.first + start] ^= 0x1;
            }
            if (pair.second != BOUNDARY_INDEX) {
                running_syndrome[pair.second + start] ^= 0x1;
            }
            visited.insert(pair.first);
            visited.insert(pair.second);
        }
        for (uint i = 0; i < res.correction.size(); i++) {
            correction[i] ^= res.correction[i];
        }
        start += params.detectors_per_round;
        end += params.detectors_per_round;
    }

    // Compare correction to observed flips.
    DecoderShotResult res = {
        params.sram_access_time / params.clock_frequency * 1e9,  // Convert to ns.
        lut_memory_overhead,
        is_logical_error(correction, syndrome, n_detectors, n_observables),
        correction
    };
    return res;
}

void
LILLIPUT::clear_stats() {
    Decoder::clear_stats();
    lut_misses = 0;
}

void
LILLIPUT::clear_lut() {
    lut.clear();
}

void
LILLIPUT::populate_lut(std::mt19937_64& rng) {
    uint n_detectors = backing_decoder_p->circuit.count_detectors();
    uint n_observables = backing_decoder_p->circuit.count_observables();

    uint n_batches = 10;
    uint32_t shots_per_batch = 10000;
    uint32_t total_shots = shots_per_batch * n_batches;

    std::array<uint32_t, 100> hamming_freq;
    hamming_freq.fill(0);
    
    uint32_t n_le = 0;  // Number of logical errors.
    for (uint b = 0; b < n_batches; b++) {
        stim::simd_bit_table sample_buffer = stim::detector_samples(
                    backing_decoder_p->circuit, shots_per_batch, false, true, rng);
        sample_buffer = sample_buffer.transposed();
#ifdef USE_OMP
#pragma omp parallel for reduction (+: n_le)
#endif
        for (uint32_t sn = 0; sn < shots_per_batch; sn++) {
            auto syndrome =
                _to_vector(sample_buffer[sn], n_detectors, n_observables);
            DecoderShotResult res = backing_decoder_p->decode_error(syndrome);
            // Update logical errors.
            n_le += res.is_logical_error ? 1 : 0;
            // Compute Hamming weight of the syndrome.
            uint hamming_weight = 0;
            for (uint i = 0; i < n_detectors; i++) {
                if (syndrome[i]) {
                    hamming_weight++;
                }
            }
#ifdef USE_OMP
#pragma omp atomic update
#endif
            hamming_freq[hamming_weight]++;
        }
    }
    // Compute Hamming weight cutoff.
    uint max_hamming_weight = 0;
    for (uint i = 1; i < hamming_freq.size(); i++) {
        if (hamming_freq[i] >= n_le) {
            max_hamming_weight = i;
        }  
    }

    uint32_t n_bits = n_detectors;
    uint32_t max_addr = (1 << n_bits) - (1 << (n_bits - max_hamming_weight));
    // For building the syndrome reference.
    uint syndrome_size = n_detectors + n_observables;
    // Initialize LUT.
    lut = std::vector<DecoderShotResult>(max_addr+1);
    is_lut_entry_valid = std::vector<uint8_t>(max_addr+1, 0);
    // Do not populate the LUT. We will do this on the fly
    // for performant simulation.
    // We only care in the number of entries we will have.
    lut_memory_overhead = max_addr * n_detectors * log2(n_detectors);
}

DecoderShotResult
LILLIPUT::get_lut_entry(const std::vector<uint8_t>& syndrome) {
    uint32_t n_detectors = backing_decoder_p->circuit.count_detectors();
    uint32_t n_observables = backing_decoder_p->circuit.count_observables();
    // TODO use the LUT to get the result.
    uint32_t addr = syndrome_to_int(syndrome, n_detectors);
    if (addr >= lut.size()) {
        lut_misses++;
        return (DecoderShotResult) {
            0.0,
            0.0,
            false,
            std::vector<uint8_t>(n_detectors, 0),
            std::map<uint,uint>()
        };
    } else {
        if (is_lut_entry_valid[addr]) {
            return lut[addr];
        } else {
            DecoderShotResult res = backing_decoder_p->decode_error(syndrome);
            lut[addr] = res;
            return res;
        }  
    }

}
