/*
 *  author: Suhas Vittal
 *  date:   9 September 2022
 * */

#include "hyperion/multi_qubit.h"

namespace qrc {
namespace hyperion {

HyperionMultiQubitSimulator::HyperionMultiQubitSimulator(
        const std::vector<stim::Circuit>& circuits,
        uint n_decoders,
        uint n_detectors_per_round,
        const HyperionParams& params)
    :
    /* Statistics */
    n_timeouts(0),
    n_overflows(0),
    n_logical_errors(0),
    max_faults(0),
    max_latency(0),
    /* Memory */
    dram(nullptr),
    memory_event_table(nullptr),
    /* Decoders */
    simulators(n_decoders),
    occupied(n_decoders, false),
    /* Metadata */
    circuits(circuits),
    graphs(circuits.size()),
    path_tables(circuits.size()),
    n_rounds(circuits[0].count_detectors() / n_detectors_per_round),
    main_clock_frequency(params.main_clock_frequency),
    dram_clock_frequency(params.dram_clock_frequency)
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
        graphs[i] = graph;
        path_tables[i] = compute_path_table(graph);
    }

    for (uint i = 0; i < n_decoders; i++) {
        uint8_t bankgroup = i % dram->config_->bankgroups;
        uint8_t bank = (i / dram->config_->bankgroups) 
                        % dram->config_->banks_per_group;
        uint32_t row_offset = ROWS_PER_QUBIT * i;

        HyperionSimulatorParams sim_params = {
            n_detectors,
            n_detectors_per_round,
            params.n_registers,
            params.bfu_fetch_width,
            params.bfu_compute_stages,
            bankgroup,
            bank,
            row_offset
        };
        simulators[i] = new HyperionSimulator(dram,
                                            memory_event_table,
                                            path_tables[0], 
                                            sim_params);
    }
}

HyperionMultiQubitSimulator::~HyperionMultiQubitSimulator() {
    dram->PrintStats();

    for (HyperionSimulator * s : simulators) {
        delete s;
    }
    delete dram;
    delete memory_event_table;
}

void
HyperionMultiQubitSimulator::benchmark(uint32_t shots, std::mt19937_64& rng) {
    const uint32_t shots_per_round = 100000;

    const fp_t min_clock_frequency = main_clock_frequency < dram_clock_frequency
                                ? main_clock_frequency : dram_clock_frequency;
    const uint32_t main_tpc = (uint32_t) 
                                (main_clock_frequency / min_clock_frequency);
    const uint32_t dram_tpc = (uint32_t)
                                (dram_clock_frequency / min_clock_frequency);

    const uint n_detectors = circuits[0].count_detectors();
    const uint n_observables = circuits[0].count_observables();

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
                service_queue.push_back(i);
            }
            if (service_queue.size() > simulators.size()) {
                n_overflows++;
            }
            if (service_queue.size() > max_faults) {
                max_faults = service_queue.size();
            }
#ifdef HMSIM_DEBUG
            std::cout << "Shot " << j << ": " 
                << service_queue.size() << " qubits have faults.\n";
#endif
            
            std::map<uint, uint> assignment_table;
            std::map<uint, std::vector<uint8_t>> syndrome_table;
            // Preload data.
            for (uint i = 0; i < simulators.size(); i++) {
                HyperionSimulator * sim = simulators[i];
                if (service_queue.empty()) {
                    sim->force_idle();
                } else {
                    uint qubit_id = service_queue.front();
                    service_queue.pop_front();
                    auto syndrome = _to_vector(sample_buffers[qubit_id][j], 
                                                n_detectors, n_observables); 
                    load_simulator(sim, syndrome, qubit_id);
                    // Mark qubit_id as being serviced by decoder i
                    assignment_table[i] = qubit_id;
                    syndrome_table[i] = syndrome;
#ifdef HMSIM_DEBUG
                    std::cout << "\tDecoder " << i << " --> qubit " 
                        << qubit_id << "\n";
#endif
                }
            }

            fp_t time_taken = 0.0; 
            uint64_t n_cycles = 0;
            uint round = 0;
            bool done;
            do {
                done = true;
                for (uint i = 0; i < dram_tpc; i++) {
                    dram->ClockTick();
                }
                for (uint i = 0; i < simulators.size(); i++) {
                    HyperionSimulator * sim = simulators[i]; 
                    if (sim->is_idle()) {
                        if (!assignment_table.count(i)
                            || assignment_table[i] == (uint)-1) 
                        {
                            continue;
                        }
                        // Check for a logical error.
                        auto matching = sim->get_matching();
                        auto correction = get_correction(matching, 
                                                        assignment_table[i]);
                        bool is_error = is_logical_error(correction, 
                                                        syndrome_table[i], 
                                                        n_detectors,
                                                        n_observables);
                        if (is_error) {
                            n_logical_errors++;
                        }
                        assignment_table[i] = (uint)-1;
                        // If the queue is empty, then all faults have
                        // been assigned.
                        if (service_queue.empty()) {
                            continue;
                        }
                        // Otherwise, we can service more events.
                        done = false;
                        uint qubit_id = service_queue.front();
                        service_queue.pop_front();
                        auto syndrome = _to_vector(sample_buffers[qubit_id][j],
                                                    n_detectors, n_observables);
                        load_simulator(sim, syndrome, qubit_id);
                        assignment_table[i] = qubit_id;
                        syndrome_table[i] = syndrome;
                        sim->sig_end_round(round);  // Need to fastforward
                                                  // the simulator to accept
                                                  // syndromes from past rounds.
#ifdef HMSIM_DEBUG
                        std::cout << "\t[ cycle = " << n_cycles 
                            << " ] Decoder " << i << " --> qubit " 
                            << qubit_id << "\n";
#endif
                    } else {
                        done = false;
                    }
                    for (uint j = 0; j < main_tpc; j++) {
                        sim->tick();
                    }
                }
                
                n_cycles++;
                time_taken = n_cycles / min_clock_frequency * 1e9;
                if (time_taken >= 1000.0) {
                    if (round < n_rounds) {
                        for (HyperionSimulator * sim : simulators) {
                            sim->sig_end_round();
                        }
#ifdef HMSIM_DEBUG
                        std::cout << "\tRound " << round << " ended.\n";
#endif
                        round++;
                        n_cycles = 0;
                        time_taken = 0.0;
                    } else {
                        done = true;
                    }
                }
            } while (!done);
            if (round < n_rounds) {
                n_cycles = 0;  // We have not reached the final round,
                               // so all our computation finishes
                               // before the final round.
            }
            // Update timing statistics.
            time_taken = n_cycles / min_clock_frequency * 1e9;
            if (time_taken >= 1000.0 && service_queue.size() > 0) {
                n_timeouts++;
                n_logical_errors += service_queue.size();
#ifdef HMSIM_DEBUG
                std::cout << "\tTimed out: " << service_queue.size()
                            << " faults left.\n";
#endif
            }

            if (time_taken > max_latency) {
                max_latency = time_taken;
            }
        }
        shots -= shots_this_round;
    }
}

void
HyperionMultiQubitSimulator::reset_stats() {
    n_timeouts = 0;
    n_overflows = 0;
    n_logical_errors = 0;
    max_faults = 0;
    max_latency = 0.0;
}

void
HyperionMultiQubitSimulator::load_simulator(
        HyperionSimulator * sim,
        const std::vector<uint8_t>& syndrome,
        uint qubit_id) 
{
    stim::Circuit circ = circuits[qubit_id];
    uint n_detectors = circ.count_detectors();
    uint n_observables = circ.count_observables();
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
#ifdef HMSIM_DEBUG
    std::cout << "\tDetectors for qubit " << qubit_id << ":";
#endif
    for (uint di = 0; di < n_detectors; di++) {
        if (syndrome[di]) {
            detector_array.push_back(di);
#ifdef HMSIM_DEBUG
            std::cout << " " << di;
#endif
        }
    }
#ifdef HMSIM_DEBUG
    std::cout << "\n";
#endif
    if (detector_array.size() & 0x1) {
        detector_array.push_back(BOUNDARY_INDEX);
    }

    sim->load_detectors(detector_array);
}

std::vector<uint8_t>
HyperionMultiQubitSimulator::get_correction(
        const std::map<uint, uint>& matching,
        uint qubit_id) 
{
    std::set<uint> visited;
    std::vector<uint8_t> correction(circuits[qubit_id].count_observables(), 0);
    auto& graph = graphs[qubit_id];
    auto& path_table = path_tables[qubit_id];

    for (auto di_dj : matching) {
        uint di = di_dj.first;
        uint dj = di_dj.second;
        if (visited.count(di) || visited.count(dj)) {
            continue;
        }
        // Check path between the two detectors.
        // This is examining the error chain.
        std::vector<uint> detector_path(path_table[di_dj].path);
        for (uint i = 1; i < detector_path.size(); i++) {
            // Get edge from decoding graph.
            auto wi = detector_path[i-1];
            auto wj = detector_path[i];
            auto edge = graph.get_edge(wi, wj);
            // The edge should exist.
            for (uint obs : edge->frames) {
                // Flip the bit.
                if (obs >= 0) {
                    correction[obs] = !correction[obs];
                }
            }
        }
        visited.insert(di);
        visited.insert(dj);
    }
    return correction;
}

}  // hyperion
}  // qrc
