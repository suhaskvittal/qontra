/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#include "gulliver.h"

namespace qrc {

Gulliver::Gulliver(const stim::Circuit circuit,
        uint n_detectors_per_round,
        const GulliverParams& params)
    :MWPMDecoder(circuit), 
    // Statistics
    n_total_accesses(0),
    n_logical_failures(0),
    max_latency(0),
    max_bfu_cycles(0),
    max_prefetch_cycles(0),
    max_hamming_weight(0),
    // Memory system
    simulator(nullptr),
    // Properties
    bfu_hw_threshold(params.bfu_hw_threshold),
    n_rounds(circuit.count_detectors()/n_detectors_per_round),
    main_clock_frequency(params.main_clock_frequency),
    // Heap initialized data (to be deleted)
    dram(nullptr),
    memory_event_table(nullptr)
{
    // Initialize Memory System.
    // Define memory event table and associated callbacks.
    memory_event_table = new std::map<addr_t, bool>();

    dramsim3::Config config(params.dram_config_file, params.log_output_directory);
    auto cb = [this](addr_t x) 
    {
        this->memory_event_table->insert_or_assign(x, true);
    };

    dram = new dramsim3::MemorySystem(params.dram_config_file, 
                                        params.log_output_directory,
                                        cb, cb);

    gulliver::GulliverSimulatorParams sim_params = {
        circuit.count_detectors()+1,
        n_detectors_per_round,
        params.n_registers,
        params.bfu_fetch_width,
        params.bfu_hw_threshold,
        0,
        0,
        0
    };
    simulator = new gulliver::GulliverSimulator(dram, 
                                    nullptr,  // Single qubit experiment doesn't
                                              // use a cache.
                                    memory_event_table,
                                    path_table, 
                                    sim_params);
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
        n_logical_failures++;
        DecoderShotResult res = {
            0.0,
            0.0,
            true,
            std::vector<uint8_t>(),
            std::map<uint, uint>()
        };
        return res;
    } else if (hw == 0) {
        DecoderShotResult res = {
            0.0,
            0.0,
            false,
            std::vector<uint8_t>(),
            std::map<uint, uint>()
        };
        return res;
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
        simulator->reset_stats();
        simulator->load_detectors(detector_array);

        uint64_t n_cycles = 0;
        uint round = 0;
        while (!simulator->is_idle()) {
            fp_t t = n_cycles / main_clock_frequency * 1e9;
            if (t > 1000 && round < n_rounds) {
                simulator->sig_end_round();
                n_cycles = 0;
                round++;
            }
            dram->ClockTick();
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
                auto wi = detector_path[i-1];
                auto wj = detector_path[i];
                auto edge = graph.get_edge(wi, wj);
                for (uint obs : edge.frames) {
                    if (obs >= 0) {
                        correction[obs] = !correction[obs];
                    }
                }
            }
            visited.insert(di);
            visited.insert(dj);
        }

        fp_t time_taken = n_cycles / main_clock_frequency * 1e9;
        if (time_taken > max_latency) {
            max_latency = time_taken;
            max_bfu_cycles = simulator->bfu_cycles;
            max_prefetch_cycles = simulator->prefetch_cycles;
            max_hamming_weight = hw;
        }

        bool is_error = 
            is_logical_error(correction, syndrome, n_detectors, n_observables);
        if (time_taken > 1000) {
            n_logical_failures++;
            is_error = true;
        }

        DecoderShotResult res = {
            time_taken,
            0.0, // TODO
            is_error,
            correction,
            matching
        };
        return res;
    }
}

}  // namespace qrc
