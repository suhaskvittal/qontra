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
    // Memory system
    cache(nullptr),
    // Properties
    n_bfu(params.n_bfu),
    n_bfu_cycles_per_add(params.n_bfu_cycles_per_add),
    bfu_hw_threshold(params.bfu_hw_threshold),
    clock_frequency(params.clock_frequency),
    _sram_cost((1 << params.cacheC))
{

    GulliverCacheParams cache_params = {
        params.cacheC,
        params.cacheS,
        params.cacheB,
        circuit.count_detectors() + 1,
        false,
        params.dram_config_file,
        params.log_output_directory
    };
    cache = new GulliverCache(cache_params);
}

Gulliver::~Gulliver() {
    delete cache;
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
    return _sram_cost;
}

uint64_t
Gulliver::dram_cost() {
    uint n_d = circuit.count_detectors();
    return n_d*n_d*sizeof(uint);
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
        // Assume initial amount of time taken is 1 cycle per
        // 1 in syndrome.
        // We will push all the active high detector ids onto
        // a "hardware" array.
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
        uint64_t n_cycles = 0;
        n_cycles += cache->prefetch(detector_array);
        uint64_t n_bfu_cycles = 0;
        std::vector<BFUResult> matchings = 
            brute_force_matchings(detector_array, n_bfu_cycles);
        // Expression:
        // ceil(sum(cycles)/#fu) + #matchings - 1
        // = (sum(cycles)-1)/#fu + 1 + #matchings - 1
        // = (sum(cycles)-1)/#fu + #matchings.
        n_bfu_cycles = ((n_bfu_cycles-1)/n_bfu) + 1;
        n_cycles += n_bfu_cycles;
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
        auto matching = best_result.matching;

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

        fp_t time_taken = (n_cycles * 1e9) / clock_frequency;  // in ns.
        if (time_taken > max_bfu_latency) {
            max_bfu_latency = time_taken;
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

std::vector<BFUResult>
Gulliver::brute_force_matchings(const std::vector<uint>& detector_array,
        uint64_t& n_cycles) 
{
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();
    uint hw = detector_array.size();
    // We try and model our hardware using the appropriate
    // data structures.
    //
    // We implement a BFU in hardware using a hardware stack.
    // and pipelined logic.
    std::stack<BFUResult> bfu_stack;
    BFUResult init = {std::map<uint,uint>(), 0.0, true};
    bfu_stack.push(init);
    // We store the results in a hardware array.
    std::vector<BFUResult> results;  // Only has to be (hw-1)! in size.

    uint res_index = 0;
    while (bfu_stack.size() > 0) {
        // Assume that popping off the hardware
        // stack is one cycle.
        BFUResult entry = bfu_stack.top();
        bfu_stack.pop();

        n_cycles++;
        // Find first unmatched detector, and match to ever other unmatched detector.
        bool found_unmatched_detector = false;
        uint first_unmatched_detector = 0;
        for (uint8_t ai = 0; ai < detector_array.size(); ai++) {
            uint di = detector_array[ai];
            // Assume searching through the matching requires
            // matching.size() cycles (one cycle per entry, assume
            // we read every single entry).
            n_cycles += entry.matching.size();
            if (entry.matching.count(di)) {
                continue;
            }

            if (found_unmatched_detector) {
                std::pair di_dj = std::make_pair(first_unmatched_detector, di);
                // Recurse with this new assignment.
                // Copy data from running result.
                std::map<uint, uint> matching(entry.matching);
                // Get data from the cache (simulated)
                n_cycles += cache->access(first_unmatched_detector, di);
                fp_t cost = entry.matching_weight
                                + path_table[di_dj].distance;
                matching[first_unmatched_detector] = di;
                matching[di] = first_unmatched_detector;
                // We perform an add operation
                // and update the matching.
                //
                // Assume the write to the matching
                // data structure is 2 cycles.
                n_cycles += n_bfu_cycles_per_add + 2;
                // Update matching.
                BFUResult new_result = {matching, cost, true};
                bfu_stack.push(new_result);
            } else {
                found_unmatched_detector = true;
                first_unmatched_detector = di;
            }
            // Assume each iteration through the detector array is one cycle.
            n_cycles++;
        } 
        // If we could not find an unmatched detector, this
        // is a perfect matching.
        if (!found_unmatched_detector) {
            results.push_back(entry);
        }
    }
    return results;
}

