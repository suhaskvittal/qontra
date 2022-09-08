/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#include "gulliver.h"

Gulliver::Gulliver(const stim::Circuit circuit,
        const GulliverParams& params)
    :MWPMDecoder(circuit), 
    // Statistics
    n_total_accesses(0),
    n_mwpm_accesses(0),
    max_bfu_latency(0),
    max_cycles(0),
    // Memory system
    simulator(nullptr),
    // Properties
    bfu_hw_threshold(params.bfu_hw_threshold),
    main_clock_frequency(params.main_clock_frequency),
    // Heap initialized data (to be deleted)
    dram(nullptr),
    memory_event_table(nullptr)
{
    // Initialize Memory System.
    // Define memory event table and associated callbacks.
    memory_event_table = new std::map<std::pair<uint, uint>, bool>();

    dramsim3::Config config(params.dram_config_file, params.log_output_directory);
    const addr_t base_address = get_base_address(0, 0, 0, &config);
    const uint n_detectors = circuit.count_detectors()+1;
    auto cb = [this, base_address, n_detectors](addr_t x) 
    {
        auto di_dj = from_address(x, base_address, n_detectors);
        this->memory_event_table->insert_or_assign(di_dj, true);
    };
    dram = new dramsim3::MemorySystem(params.dram_config_file, 
                                        params.log_output_directory,
                                        cb, cb);
    std::map<std::pair<uint, uint>, fp_t> weight_table;
    for (auto res : path_table) {
        auto di_dj = res.first;
        DijkstraResult r = res.second;
        weight_table[di_dj] = r.distance;
    }

    GulliverSimulatorParams sim_params = {
        n_detectors,
        params.n_registers,
        params.bfu_fetch_width,
        params.bfu_hw_threshold,
        0,
        0,
        0
    };
    simulator = new GulliverSimulator(dram, this->memory_event_table,
                                    weight_table, sim_params);
                    
}

Gulliver::~Gulliver() {
    dram->PrintStats();

    delete simulator;
    delete dram;
    delete memory_event_table;
}

std::string
Gulliver::name() {
    return "Gulliver";
}

bool
Gulliver::is_software() {
    return false;
}

uint64_t
Gulliver::sram_cost() {
    return 0;   // TODO
}

uint64_t
Gulliver::dram_cost() {
    uint n_d = circuit.count_detectors();
    return n_d*(n_d-1)*sizeof(uint)/2;  // In bytes.
}

DecoderShotResult
Gulliver::decode_error(const std::vector<uint8_t>& syndrome) {
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();
    // Compute Hamming weight.
    // Don't count the observables.
    uint hw = std::accumulate(syndrome.begin(), syndrome.end()-n_observables, 0);
    n_total_accesses++;
    if (hw > bfu_hw_threshold) {
        n_mwpm_accesses++;
        return MWPMDecoder::decode_error(syndrome);
    } else {
        // Invoke BFUs.
        // We use a simulator to count the number of cycles.
        std::vector<uint> detector_array;
        for (uint di = 0; di < n_detectors; di++) {
            if (syndrome[di]) {
                detector_array.push_back(di);
            }
        }
        if (detector_array.size() & 0x1) {
            // Add boundary.
            detector_array.push_back(BOUNDARY_INDEX);
        }
        // Run simulator.
        simulator->load(detector_array);
        uint64_t n_cycles = 0;
        while (!simulator->is_idle()) {
            simulator->tick();
            n_cycles++;
        }
        // Get matching from simulator.
        auto matching = simulator->get_matching();

        std::vector<uint8_t> correction(n_observables, 0);
        std::set<uint> visited;
        for (auto assignment : matching) {
            uint di = assignment.first;
            uint dj = assignment.second;
            if (visited.count(di) || visited.count(dj)) {
                continue;
            }
            // For explanation, see MWPMDecoder function.
            std::pair<uint, uint> di_dj = std::make_pair(di, dj);
            std::vector<uint> detector_path(path_table[di_dj].path);
            for (uint i = 1; i < detector_path.size(); i++) {
                // Get edge from decoding graph.
                auto wi = graph.get(detector_path[i-1]);
                auto wj = graph.get(detector_path[i]);
                auto edge = boost::edge(wi, wj, graph.base);
                for (uint obs : graph.base[edge.first].frames) {
                    if (obs >= 0) {
                        correction[obs] = !correction[obs];
                    }
                }
            }
            visited.insert(di);
            visited.insert(dj);
        }

        fp_t time_taken = n_cycles / main_clock_frequency * 1e9;
        if (time_taken > max_bfu_latency) {
            max_bfu_latency = time_taken;
            max_cycles = n_cycles;
        }
        DecoderShotResult res = {
            time_taken,
            0.0, // TODO
            is_logical_error(correction, syndrome, n_detectors, n_observables),
            correction,
            matching
        };
        return res;
    }
}

/* FOR FUTURE REFERENCE, DO NOT DELETE
 *
 *
 *
std::map<uint, uint>
Gulliver::predecode(const std::vector<uint>& detector_array, 
        GulliverCycles& n_cycles) 
{
    uint64_t n_proc_cycles = 0;
    GulliverCycles n_mem_cycles;  // sram and dram accesses.
    
    typedef std::pair<uint, fp_t> PairingEntry;
    std::map<uint, PairingEntry> best_pairing; 
    std::map<uint, uint> n_pairings;
    for (uint i = 0; i < detector_array.size(); i++) {
        uint di = detector_array[i];
        if (n_pairings.count(di) == 0) {
            n_pairings[di] = 0;
        }
        auto vi = graph.get(di);
        // Compute the min neighbor that is in the detector array.
        for (uint j = i+1; j < detector_array.size(); j++) {
            uint dj = detector_array[j];
            if (n_pairings.count(dj) == 0) {
                n_pairings[dj] = 0;
            }
            auto vj = graph.get(dj);
            // Check if di and dj are connected.
            // Assume this is an access to some SRAM data structure.
            //
            // This can be a global data structure across all logical qubits.
            // Or even just a logic circuit.
            auto edge = boost::edge(vi, vj, graph.base);
            n_proc_cycles++;
            if (!edge.second) {
                continue;
            }
            // If they are adjacent, then access the weight from memory.
            fp_t w = graph.base[edge.first].edge_weight;
            n_mem_cycles += memsys->access(di, dj, false);

            if (best_pairing.count(di)) {
                PairingEntry prev = best_pairing[di];
                if (prev.second > w) {
                    best_pairing[di] = std::make_pair(dj, w);
                    n_pairings[dj]++;
                    n_pairings[prev.first]--;
                }
            } else {
                best_pairing[di] = std::make_pair(dj, w);
                n_pairings[dj]++;
            }

            if (best_pairing.count(dj)) {
                PairingEntry prev = best_pairing[dj];
                if (prev.second > w) {
                    best_pairing[dj] = std::make_pair(di, w);
                    n_pairings[di]++;
                    n_pairings[prev.first]--;
                }
            } else {
                best_pairing[dj] = std::make_pair(di, w);
                n_pairings[di]++;
            }
            // Assume one cycle for comparisons.
            n_proc_cycles++;
        }
    }
    // Consolidate the matching.
    std::map<uint, uint> matching;
    for (auto kv_entry : best_pairing) {
        uint di = kv_entry.first;
        uint dj = kv_entry.second.first;
        if (matching.count(di)) {
            continue;
        }
        if (n_pairings[di] > 1 || n_pairings[dj] > 1) {
            continue;
        }
        if (di == best_pairing[dj].first) {
            matching[di] = dj;
            matching[dj] = di;
        }
        n_proc_cycles++;
    }
    n_cycles.onchip += n_proc_cycles;
    n_cycles += n_mem_cycles;
    return matching;
}
*/
