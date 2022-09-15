/*
 *  author: Suhas Vittal
 *  date:   9 September 2022
 * */

#include "gulliver/multi_qubit.h"

namespace qrc {
namespace gulliver {

GulliverMultiQubitSimulator::GulliverMultiQubitSimulator(
        const std::vector<stim::Circuit>& circuits,
        uint n_decoders,
        uint n_detectors_per_round,
        const GulliverParams& params)
    :
    /* Statistics */
    n_timeouts(0),
    n_overflows(0),
    n_uncomputable(0),
    max_latency(0),
    /* Memory */
    dram(nullptr),
    memory_event_table(nullptr),
    /* Decoders */
    simulators(n_decoders),
    occupied(n_decoders, false),
    /* Metadata */
    circuits(circuits),
    path_tables(circuits.size()),
    main_clock_frequency(params.main_clock_frequency)
{
    // Initialize memory first.
    memory_event_table = new std::map<addr_t, bool>(); 
    // Note that every circuit, provided that they are all
    // the same code and code distance, have the same
    // adjacency matrix.
    auto cb = [this](addr_t x) 
    {
        this->memory_event_table->insert_or_assign(x, true);
    };
    dram = new dramsim3::MemorySystem(params.dram_config_file,
                                        params.log_output_directory,
                                        cb, cb);
    for (uint i = 0; i < circuits.size(); i++) {
        DecodingGraph graph = to_decoding_graph(circuits[i]);
        path_tables[i] = compute_path_table(graph);
    }

    for (uint i = 0; i < n_decoders; i++) {
        uint8_t bankgroup = i % dram->config_->bankgroups;
        uint8_t bank = (i / dram->config_->bankgroups) 
                        % dram->config_->banks_per_group;
        uint32_t row_offset = ROWS_PER_QUBIT * i;

        GulliverSimulatorParams sim_params = {
            circuits[0].count_detectors()+1,
            n_detectors_per_round,
            params.n_registers,
            params.bfu_fetch_width,
            params.bfu_hw_threshold,
            bankgroup,
            bank,
            row_offset
        };
        simulators[i] = new GulliverSimulator(dram,
                                            memory_event_table,
                                            path_tables[0], 
                                            sim_params);
    }
}

GulliverMultiQubitSimulator::~GulliverMultiQubitSimulator() {
    dram->PrintStats();

    for (GulliverSimulator * s : simulators) {
        delete s;
    }
    delete dram;
    delete memory_event_table;
}

void
GulliverMultiQubitSimulator::benchmark(uint32_t shots, std::mt19937_64& rng) {
    const uint32_t shots_per_round = 100000;

    while (shots > 0) {
        uint32_t shots_this_round = 
            shots > shots_per_round ? shots_per_round : shots;

        std::vector<stim::simd_bit_table> sample_buffers;
        for (uint i = 0; i < circuits.size(); i++) {
            auto buf = stim::detector_samples(circuits[i], shots_this_round,
                                                false, true, rng);
            buf = buf.transposed();
            sample_buffers.push_back(buf);
        }

        for (uint32_t j = 0; j < shots_this_round; j++) {
            uint next_available_decoder = 0;
            for (uint i = 0; i < circuits.size(); i++) {
                if (!sample_buffers[i][j].not_zero()) {
                    continue;
                }
                if (next_available_decoder >= simulators.size()) {
                    next_available_decoder = 0;
                    n_overflows++;
                    break;
                }

                stim::Circuit circ = circuits[i];
                uint n_detectors = circ.count_detectors();
                uint n_observables = circ.count_observables();
                auto syndrome = 
                    _to_vector(sample_buffers[i][j], n_detectors, n_observables);
                // Load data into simulator.
                auto path_table = path_tables[i];
                uint8_t bankgroup = i % dram->config_->bankgroups;
                uint8_t bank = (i / dram->config_->bankgroups) 
                                % dram->config_->banks_per_group;
                uint32_t row_offset = ROWS_PER_QUBIT * i;
                GulliverSimulator * sim = simulators[next_available_decoder];
                sim->load_path_table(path_table);
                sim->load_base_address(bankgroup, bank, row_offset);
                // Finally, load the problem into the simulator.
                std::vector<uint> detector_array;
                for (uint di = 0; di < n_detectors; di++) {
                    if (syndrome[di]) {
                        detector_array.push_back(di);
                    }
                }

                if (!sim->load_detectors(detector_array)) {
                    n_uncomputable++;
                }
                next_available_decoder++;
            }

            uint64_t n_cycles = 0;
            bool done;
            do {
                done = true;
                for (uint i = 0; i < next_available_decoder; i++) {
                    GulliverSimulator * sim = simulators[i];
                    if (!sim->is_idle()) {
                        done = false;
                        sim->tick();
                    }
                }
                n_cycles++;
            } while (!done);
            fp_t time_taken = n_cycles / main_clock_frequency * 1e9;
            if (time_taken > 1000.0) {
                n_timeouts++;
            }

            if (time_taken > max_latency) {
                max_latency = time_taken;
            }
        }
        shots -= shots_this_round;
    }
}

void
GulliverMultiQubitSimulator::reset_stats() {
    n_timeouts = 0;
    n_overflows = 0;
    n_uncomputable = 0;
    max_latency = 0.0;
}

}  // gulliver
}  // qrc
