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
    max_cycles(),
    // Memory system
    memsys(nullptr),
    // Properties
    n_bfu(params.n_bfu),
    n_bfu_cycles_per_add(params.n_bfu_cycles_per_add),
    bfu_hw_threshold(params.bfu_hw_threshold),
    main_clock_frequency(params.main_clock_frequency),
    dram_clock_frequency(params.dram_clock_frequency)
{

    GulliverMemoryParams mem_params = {
        params.n_sram_table_entries,
        circuit.count_detectors() + 1,
        params.dram,
        params.bankgroup,
        params.bank,
        params.row_offset,
        params.dram_config_file,
        params.log_output_directory
    };
    memsys = new GulliverMemory(mem_params);
}

Gulliver::~Gulliver() {
    delete memsys;
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
        // Invoke BFUs. As this is hardware, we don't use system time
        // to track time taken.
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
        GulliverCycles n_cycles;
        n_cycles += memsys->prefetch(detector_array);
        // First do some predecoding.
        auto matching = predecode(detector_array, n_cycles);
        // Remove entries from the detector array that have been
        // matched.
        for (auto it = detector_array.begin(); it != detector_array.end(); ) {
            if (matching.count(*it)) {
                it = detector_array.erase(it);
            } else {
                it++;
            }
        }
        std::vector<BFUResult> matchings = 
            brute_force_matchings(detector_array, n_cycles);
        // Choose the best one -- we assume this takes
        // #matchings-1 cycles though this could be done in 1
        // cycle with enough comparator gates.
        BFUResult best_result = matchings[0];
        for (uint i = 1; i < matchings.size(); i++) {
            auto m = matchings[i];
            if (m.matching_weight < best_result.matching_weight) {
                best_result = m;
            }
        }
        // Compute logical correction.
        for (auto kv_entry : best_result.matching) {
            matching[kv_entry.first] = kv_entry.second; 
        }

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

        fp_t time_taken = 
            ((n_cycles.onchip / main_clock_frequency)
             + (n_cycles.dram / dram_clock_frequency)) * 1e9;  // in ns.
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

std::vector<BFUResult>
Gulliver::brute_force_matchings(const std::vector<uint>& detector_array,
        GulliverCycles& n_cycles) 
{
    uint64_t n_proc_cycles = 0;
    GulliverCycles n_mem_cycles;  // sram and dram.

    uint hw = detector_array.size();
    // We try and model our hardware using the appropriate
    // data structures.
    //
    // We implement a BFU in hardware using a hardware stack.
    // and pipelined logic.
    uint max_stack_size = 0;  // We track this to see number of BFUs being used.
    std::stack<BFUResult> bfu_stack;
    BFUResult init = {std::map<uint,uint>(), 0.0, true};
    bfu_stack.push(init);
    // We store the results in a hardware array.
    std::vector<BFUResult> results;  // Only has to be (hw-1)! in size.

    uint res_index = 0;
    while (bfu_stack.size() > 0) {
        if (bfu_stack.size() > max_stack_size) {
            max_stack_size = bfu_stack.size();
        }
        // Assume that popping off the hardware
        // stack is one cycle.
        BFUResult entry = bfu_stack.top();
        bfu_stack.pop();
        n_proc_cycles++;
        // Find first unmatched detector, and match to ever other unmatched detector.
        bool found_unmatched_detector = false;
        uint first_unmatched_detector = 0;
        for (uint8_t ai = 0; ai < detector_array.size(); ai++) {
            uint di = detector_array[ai];
            n_proc_cycles++;
            if (entry.matching.count(di)) {
                continue;
            }

            if (found_unmatched_detector) {
                std::pair di_dj = std::make_pair(first_unmatched_detector, di);
                // Recurse with this new assignment.
                // Copy data from running result.
                std::map<uint, uint> matching(entry.matching);
                // Get data from the memory hierarchy
                n_mem_cycles += memsys->access(first_unmatched_detector, di, true);
                fp_t cost = entry.matching_weight
                                + path_table[di_dj].distance;
                matching[first_unmatched_detector] = di;
                matching[di] = first_unmatched_detector;
                // We perform an add operation
                // and update the matching.
                //
                // Assume the write to the matching
                // data structure is 2 cycles.
                n_proc_cycles += n_bfu_cycles_per_add + 2;
                // Update matching.
                BFUResult new_result = {matching, cost, true};
                bfu_stack.push(new_result);
            } else {
                found_unmatched_detector = true;
                first_unmatched_detector = di;
            }
        } 
        // If we could not find an unmatched detector, this
        // is a perfect matching.
        if (!found_unmatched_detector) {
            results.push_back(entry);
        }
    }
    n_cycles.onchip += n_proc_cycles;
    n_cycles += n_mem_cycles;
    // Observe that accesses to the registers can be concurrent,
    // so we can assume that any step, all BFUs are using the
    // register file.
    uint n_bfu_used = n_bfu < max_stack_size ? n_bfu : max_stack_size;
    n_cycles.onchip = (n_cycles.onchip-1)/n_bfu_used + 1;
    return results;
}

