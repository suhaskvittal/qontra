/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#include "quarch_gulliver.h"

GulliverDecoder::GulliverDecoder(const stim::Circuit circuit,
        const GulliverParams& params)
    :MWPMDecoder(circuit), 
    n_bfu(params.n_bfu),
    n_bfu_cycles_per_add(params.n_bfu_cycles_per_add),
    bfu_hw_threshold(params.bfu_hw_threshold),
    clock_frequency(params.clock_frequency)
{}

std::string
GulliverDecoder::name() {
    return "Gulliver";
}

bool
GulliverDecoder::is_software() {
    return false;
}

DecoderShotResult
GulliverDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();
    // Compute Hamming weight.
    // Don't count the observables.
    uint hw = std::accumulate(syndrome.begin(), syndrome.end()-n_observables, 0);
    if (hw > bfu_hw_threshold) {
        return MWPMDecoder::decode_error(syndrome);
    } else {
        // Invoke BFUs. As this is hardware, we don't use system time
        // to track time taken.
        // Assume initial amount of time taken is 1 cycle per
        // bit in syndrome.
        GulliverDecoder::BFUResult init = 
            std::make_tuple(std::map<uint,uint>(), 0, n_detectors);
        std::vector<GulliverDecoder::BFUResult> matchings = 
            brute_force_matchings(syndrome, init);
        // Choose the best one -- we assume this takes
        // #matchings-1 cycles though this could be down in 1
        // cycle with enough comparator gates.
        GulliverDecoder::BFUResult best_result = matchings[0];
        for (uint i = 1; i < matchings.size(); i++) {
            auto m = matchings[i];
            if (std::get<1>(m) < std::get<1>(best_result)) {
                best_result = m;
            }
        }
        // Compute logical correction.
        auto matching = std::get<0>(best_result);
        // Expression:
        // ceil(sum(cycles)/#fu) + #matchings - 1
        // = (sum(cycles)-1)/#fu + 1 + #matchings - 1
        // = (sum(cycles)-1)/#fu + #matchings.
        uint32_t n_cycles = 
            ((std::get<2>(best_result)-1)/n_bfu_cycles_per_add) + matchings.size();

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

std::vector<GulliverDecoder::BFUResult>
GulliverDecoder::brute_force_matchings(const std::vector<uint8_t>& syndrome, 
        const GulliverDecoder::BFUResult& running_result)
{
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();
    uint hw = std::accumulate(syndrome.begin(), syndrome.end()-n_observables, 0);
    // Find first unmatched detector, and match to ever other unmatched detector.
    bool found_unmatched_detector = false;
    uint first_unmatched_detector = 0;
    std::vector<GulliverDecoder::BFUResult> results;
    for (uint di = 0; di < n_detectors; di++) {
        if (!syndrome[di] || std::get<0>(running_result).count(di)) {
            continue;
        }

        if (found_unmatched_detector) {
            std::pair di_dj = std::make_pair(first_unmatched_detector, di);
            // Recurse with this new assignment.
            // Copy data from running result.
            std::map<uint, uint> matching(std::get<0>(running_result));
            fp_t cost = std::get<1>(running_result)
                            + path_table[di_dj].distance;
            // Assume it took hw - 2*#matches to get to this point, plus
            // an addition operation. This code is unoptimized.
            uint32_t n_cycles = std::get<2>(running_result) 
                            + (hw - matching.size())
                            + n_bfu_cycles_per_add
                            + (hw & 0x1 ? 1 : 0);  // Add 1 if hw is odd (boundary)
            // Update matching.
            matching[first_unmatched_detector] = di;
            matching[di] = first_unmatched_detector;
            const GulliverDecoder::BFUResult& new_result
                = std::make_tuple(matching, cost, n_cycles); 
            auto sub_results = brute_force_matchings(syndrome, new_result);
            for (auto r : sub_results) {
                results.push_back(r); 
            }
        } else {
            found_unmatched_detector = true;
            first_unmatched_detector = di;
        }
    }
    if (found_unmatched_detector) {
        if (hw & 0x1) {
            // Also try matching to the boundary. 
            std::pair di_dj = 
                std::make_pair(first_unmatched_detector, BOUNDARY_INDEX);
            std::map<uint, uint> matching(std::get<0>(running_result));
            fp_t cost = std::get<1>(running_result)
                            + path_table[di_dj].distance;
            uint32_t n_cycles = std::get<2>(running_result) 
                            + (hw - matching.size())
                            + n_bfu_cycles_per_add
                            + (hw & 0x1 ? 1 : 0);
            matching[first_unmatched_detector] = BOUNDARY_INDEX;
            matching[BOUNDARY_INDEX] = first_unmatched_detector;
            const GulliverDecoder::BFUResult& new_result
                = std::make_tuple(matching, cost, n_cycles); 
            auto sub_results = brute_force_matchings(syndrome, new_result);
            for (auto r : sub_results) {
                results.push_back(r); 
            }
        }
        return results;
    } else {
        // The running result is a complete matching.
        return std::vector<GulliverDecoder::BFUResult>{running_result};
    }
}

