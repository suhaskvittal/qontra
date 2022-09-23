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
    caches(),
    memory_event_table(nullptr),
    /* Decoders */
    simulators(n_decoders),
    occupied(n_decoders, false),
    /* Metadata */
    circuits(circuits),
    path_tables(circuits.size()),
    n_rounds(circuits[0].count_detectors() / n_detectors_per_round),
    main_clock_frequency(params.main_clock_frequency)
{
    uint n_detectors = circuits[0].count_detectors()+1;
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

        QubitCache * cache = new QubitCache(params.n_cache_supertags,
                                params.n_cache_sets,
                                params.n_cache_lines,
                                n_detectors - n_detectors_per_round);
        caches.push_back(cache);

        GulliverSimulatorParams sim_params = {
            n_detectors,
            n_detectors_per_round,
            params.n_registers,
            params.bfu_fetch_width,
            params.bfu_hw_threshold,
            bankgroup,
            bank,
            row_offset
        };
        simulators[i] = new GulliverSimulator(dram,
                                            nullptr,
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
    for (QubitCache * c : caches) {
        delete c;
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
            std::deque<uint> service_queue;

            for (uint i = 0; i < circuits.size(); i++) {
                if (!sample_buffers[i][j].not_zero()) {
                    continue;
                }
                if (service_queue.size() >= simulators.size()) {
                    n_overflows++;
                }
                service_queue.push_back(i);
            }
#ifdef GSIM_DEBUG
            std::cout << "Shot " << j << ": " 
                << service_queue.size() << " qubits have faults.\n";
#endif
            
            // Preload data.
            for (uint i = 0; i < simulators.size(); i++) {
                GulliverSimulator * sim = simulators[i];
                if (service_queue.empty()) {
                    sim->force_idle();
                } else {
                    // Compute data
                    uint qubit_id = service_queue.front();
                    service_queue.pop_front();
                    load_simulator(sim, sample_buffers[qubit_id][j], qubit_id);
                }
            }

            fp_t time_taken = 0.0; 
            uint64_t n_cycles = 0;
            uint round = 0;
            bool done;
            do {
#ifdef GSIM_DEBUG
                if (n_cycles == 1001) {
                    std::cout << "\tTimeout: " 
                        << service_queue.size() << " qubits unserviced.\n";
                }
#endif
                done = true;
                dram->ClockTick();
                for (uint i = 0; i < simulators.size(); i++) {
                    GulliverSimulator * sim = simulators[i]; 
                    if (sim->is_idle() && !service_queue.empty()) {
                        done = false;

                        uint qubit_id = service_queue.front();
                        service_queue.pop_front();
                        load_simulator(sim, sample_buffers[qubit_id][j], qubit_id);
                        sim->sig_end_round(round);
#ifdef GSIM_DEBUG
                        std::cout << "Simulator " << i << 
                            " finished early, loading new syndrome.\n";
#endif
                    } else if (!sim->is_idle()) {
                        done = false;
                    } else {
                        continue;
                    }
#ifdef GSIM_DEBUG
                    if (n_cycles == 1001) {
                        std::cout << "Simulator " << i << " is not finished\n";
                    }
#endif
                    sim->tick();
                }
                
                n_cycles++;
                time_taken = n_cycles / main_clock_frequency * 1e9;
                if (time_taken > 1000.0 && round < n_rounds) {
                    for (GulliverSimulator * sim : simulators) {
                        sim->sig_end_round();
                    }
                    round++;
                    n_cycles = 0;
                    time_taken = 0.0;
                }
            } while (!done);
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

bool
GulliverMultiQubitSimulator::load_simulator(
        GulliverSimulator * sim,
        const stim::simd_bits_range_ref& event,
        uint qubit_id) 
{
    stim::Circuit circ = circuits[qubit_id];
    uint n_detectors = circ.count_detectors();
    uint n_observables = circ.count_observables();
    auto syndrome = _to_vector(event, n_detectors, n_observables);
    // Load data into simulator.
    auto path_table = path_tables[qubit_id];
    uint8_t bankgroup = qubit_id % dram->config_->bankgroups;
    uint8_t bank = (qubit_id / dram->config_->bankgroups)
                    % dram->config_->banks_per_group;
    uint32_t row_offset = ROWS_PER_QUBIT * qubit_id;

    sim->load_qubit_number(qubit_id);
    sim->load_path_table(path_table);
    sim->load_base_address(bankgroup, bank, row_offset);
    // Load the problem into the simulator.
    std::vector<uint> detector_array;
    for (uint di = 0; di < syndrome.size(); di++) {
        if (syndrome[di]) {
            detector_array.push_back(di);
        }
    }

    bool computable = sim->load_detectors(detector_array);
    if (!computable) {
        n_uncomputable++;
    }
    return computable;
}

}  // gulliver
}  // qrc
