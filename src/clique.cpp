/*
 *  author: Suhas Vittal
 *  date:   25 August 2022
 * */

#include "clique.h"

CliqueDecoder::CliqueDecoder(const stim::Circuit& circ, const CliqueParams& params)
    :MWPMDecoder(circ),
    n_mwpm_accesses(0),
    n_total_accesses(0),
    n_cycles_AND(params.n_cycles_AND),
    n_cycles_NOT(params.n_cycles_NOT),
    n_cycles_XOR(params.n_cycles_XOR),
    clock_frequency(params.clock_frequency),
    detectors_per_round(params.detectors_per_round),
    first_round_degree_table()
{
    // Initialize first round degree table.
    boost::graph_traits<decoding_graph_base>::vertex_iterator vi_p, vf_p;
    boost::tie(vi_p, vf_p) = boost::vertices(graph.base);
    for (; vi_p != vf_p; vi_p++) {
        auto vi = *vi_p;
        auto di = graph.base[vi].detector;
        if (di >= detectors_per_round) {
            continue; 
        }
        uint deg = 0;

        boost::graph_traits<decoding_graph_base>::adjacency_iterator ai_p, af_p;
        boost::tie(ai_p, af_p) = boost::adjacent_vertices(vi, graph.base);
        for (; ai_p != af_p; ai_p++) {
            auto vj = *ai_p;
            auto dj = graph.base[vj].detector;
            if (dj < detectors_per_round) {
                deg++; 
            }
        }
        first_round_degree_table[di] = deg;
    }
}

std::string
CliqueDecoder::name() {
    return "CliqueDecoder";
}

bool
CliqueDecoder::is_software() {
    return false;
}

DecoderShotResult
CliqueDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();

    uint64_t n_cycles_for_clique = 0;

    bool is_zero = true;
    // Read syndrome and compute which bits are active high.
    std::vector<uint8_t> parity_bits(detectors_per_round);
    uint di_mod_dpr = 0;
    for (uint di = 0; di < n_detectors; di++) {
        is_zero &= (syndrome[di] == 0);
        if (di < detectors_per_round) {
            // Count number of detectors adjacent to this one in
            // the first round.
            uint deg = first_round_degree_table[di];
            parity_bits[di] = syndrome[di] ^ (deg & 0x1);
        } else {
            // This is a measurement from a future round.
            uint8_t syndrome_stayed_extra_cycle
                = syndrome[di] ^ syndrome[di - detectors_per_round];
            parity_bits[di_mod_dpr] &= !syndrome_stayed_extra_cycle;
            // Update cycles.
            n_cycles_for_clique += n_cycles_XOR + n_cycles_NOT + n_cycles_AND;
        }
        // Update di_mod_dpr.
        // It goes back to 0 if it reaches detectors_per_round (dpr).
        if (++di_mod_dpr == detectors_per_round) {
            di_mod_dpr = 0;
        }
    }

    if (is_zero) {
        return (DecoderShotResult) {
            n_cycles_for_clique * 1e9 / clock_frequency,
            0.0,
            false,
            std::vector<uint8_t>(n_observables, 0)
        };
    }
    n_total_accesses++;
    std::vector<uint8_t> correction(n_observables, 0);
    std::vector<uint8_t> post_clique_syndrome(syndrome);
    // Detect if we need to use MWPM.
    bool goto_complex = false;
    std::vector<uint8_t> neighbor_parity_bits(detectors_per_round, 0);
    for (uint di = 0; di < detectors_per_round; di++) {
        auto vi = graph.get(di);
        // Examine adjacency list and see if clique condition
        // is maintained.
        //
        // If di && !XOR(all neighbors) then GOTO COMPLEX.
        // else, correct error.
        boost::graph_traits<decoding_graph_base>::adjacency_iterator ai_p, af_p;
        boost::tie(ai_p, af_p) = boost::adjacent_vertices(vi, graph.base);
        uint error_dj = di;  // If not a complex error,
                              // then this will contain the faulting
                              // detector.
        for (; ai_p != af_p; ai_p++) {
            auto dj = graph.base[*ai_p].detector;
            if (dj >= detectors_per_round) {
                continue;  // This one doesn't count -- it is not in the first round.
            }
            neighbor_parity_bits[di] ^= parity_bits[dj];
            // Update cycles.
            n_cycles_for_clique += n_cycles_XOR;
            if (parity_bits[di] && parity_bits[dj]) {
                error_dj = dj;
            }
        }

        bool is_complex = parity_bits[di] && !neighbor_parity_bits[di];
        n_cycles_for_clique += n_cycles_AND + n_cycles_NOT;
        goto_complex |= is_complex;
        if (!is_complex && parity_bits[di]) {
            // Go and fix error.
            auto vj = graph.get(error_dj);
            auto edge = boost::edge(vi, vj, graph.base);
            for (uint obs : graph.base[edge.first].frames) {
                if (obs >= 0) {
                    correction[obs] = !correction[obs];
                }
            }
            // Go update the syndrome as well.
            // Clear all bits in all rounds corresponding to the
            // faults.
            for (uint ddi = di; ddi < n_detectors; ddi += detectors_per_round) {
                post_clique_syndrome[ddi] = 0;
            }
            for (uint ddj = error_dj; ddj < n_detectors; 
                    ddj += detectors_per_round)
            {
                post_clique_syndrome[ddj] = 0;
            }
            // Reset parity bits as well.
            parity_bits[di] = first_round_degree_table[di] & 0x1;
            parity_bits[error_dj] = first_round_degree_table[error_dj] & 0x1;
        }
    }
    
    DecoderShotResult res = {
        n_cycles_for_clique * 1e9 / clock_frequency,
        0.0,
        is_logical_error(correction, syndrome, n_detectors, n_observables),
        correction
    };
    if (goto_complex) {
        DecoderShotResult subres = MWPMDecoder::decode_error(post_clique_syndrome); 
        res.execution_time += subres.execution_time;
        res.memory_overhead += subres.memory_overhead;
        for (uint i = 0; i < n_observables; i++) {
            res.correction[i] ^= subres.correction[i];
        }
        res.is_logical_error = is_logical_error(res.correction,
                                syndrome, n_detectors, n_observables);
        n_mwpm_accesses++;
    }
    return res;
}
