/*
 *  author: Suhas Vittal
 *  date:   2 October 2022
 * */

#include "gulliver/bdc_decoder.h"
#include <limits>

namespace qrc {
namespace gulliver {

BDCDecoder::BDCDecoder(const stim::Circuit& circ)
    :MWPMDecoder(circ)
{}

DecoderShotResult
BDCDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
    // Log start time.
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();
#ifdef __APPLE__
    auto start_time = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
#else
    struct timespec start_time_data;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start_time_data);
    auto start_time = start_time_data.tv_nsec;
#endif
    bool syndrome_is_even = true;
    std::vector<uint> detector_list;
    for (uint di = 0; di < n_detectors; di++) {
        auto syndrome_bit = syndrome[di];
        if (syndrome_bit) {
            syndrome_is_even = !syndrome_is_even;
            detector_list.push_back(di);
        }
    }
    if (!syndrome_is_even) {
        detector_list.push_back(BOUNDARY_INDEX);
    }
#ifdef BDC_DEBUG
    std::cout << "===========================================\n";
#endif
    auto matching = bdc_mwpm(detector_list);

    auto correction = get_correction_from_matching(matching);
    bool is_error = is_logical_error(correction, syndrome, 
            n_detectors, n_observables);

#ifdef BDC_DEBUG
    auto correct_result = MWPMDecoder::decode_error(syndrome);
    auto mwpm_matching = correct_result.matching;
    auto mwpm_correction = correct_result.correction;
    bool mwpm_is_error = is_logical_error(mwpm_correction,
            syndrome, n_detectors, n_observables);

    if (is_error && !mwpm_is_error) {
        // Print out weights.
        for (uint i = 0; i < detector_list.size(); i++) {
            uint di = detector_list[i];
            for (uint j = i+1; j < detector_list.size(); j++) {
                uint dj = detector_list[j];
                fp_t w = path_table[std::make_pair(di, dj)].distance;
                std::cout << "Weight(" << di << "," << dj << ") = " << w << "\n";
            }
        }
        // Check against MWPM.
        fp_t bdc_weight = 0.0;
        fp_t mwpm_weight = 0.0;
        std::cout << "BDC matching:\n";
        for (auto kv_pair : matching) {
            auto w = path_table[kv_pair].distance * 0.5;
            bdc_weight += w;
            std::cout << "\t" << kv_pair.first << " --> " << kv_pair.second << "\n";
            std::cout << "\t\tweight = " <<  w << "\n";
        }
        std::cout << "MWPM matching:\n";
        for (auto kv_pair : mwpm_matching) {
            auto w = path_table[kv_pair].distance * 0.5;
            mwpm_weight += w;
            std::cout << "\t" << kv_pair.first << " --> " << kv_pair.second << "\n";
            std::cout << "\t\tweight = " <<  w << "\n";
        }
        std::cout << "hamming weight: " << detector_list.size()
                << ", weights: " << bdc_weight << "," << mwpm_weight << "\n";
    }
#endif

#ifdef __APPLE__
    auto end_time = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
#else
    struct timespec end_time_data;
    clock_gettime(CLOCK_MONOTONIC_RAW, &end_time_data);
    auto end_time = end_time_data.tv_nsec;
#endif
    auto time_taken = end_time-start_time;
    DecoderShotResult res = {
        time_taken,
        0.0, // TODO
        is_error,
        correction,
        matching
    };
    return res;
}

std::map<uint, uint>
BDCDecoder::bdc_mwpm(const std::vector<uint>& detector_array) {
    wgt_t min_weight = std::numeric_limits<wgt_t>::max();
    // First compute the maximum weight.
    for (uint i = 0; i < detector_array.size(); i++) {
        uint di = detector_array[i];
        for (uint j = i+1; j < detector_array.size(); j++) {
            uint dj = detector_array[j];
            fp_t raw_weight = path_table[std::make_pair(di, dj)].distance;
            wgt_t w = (wgt_t) (raw_weight * MWPM_INTEGER_SCALE);

            if (w < min_weight) {
                min_weight = w;
            }
        }
    }

    typedef std::pair<uint, bool>           vertex_t;
    typedef std::pair<vertex_t, vertex_t>   edge_t;

    std::map<vertex_t, vertex_t> bp_matching;
    std::map<vertex_t, wgt_t> u_table;
    for (uint d : detector_array) {
        vertex_t v1 = std::make_pair(d, true);
        vertex_t v2 = std::make_pair(d, false);
        u_table[v1] = min_weight;
        u_table[v2] = 0;
    } 
    
    bool new_stage = true;
    std::map<vertex_t, uint8_t> label_table;  
    std::map<vertex_t, vertex_t> prev;
    while (bp_matching.size() < 2*detector_array.size()) {
        // A new stage only begins if a matching had occurred successfully.
        // At the start of every stage, we reset all data structures.
        std::set<vertex_t> visited;
        std::deque<vertex_t> vertices;
        if (new_stage) {
#ifdef BDC_DEBUG
            std::cout << "STAGE START: matching has size "
                << bp_matching.size() << "\n";
#endif
            label_table.clear();
            visited.clear();
            // Make initial label for vertex and find initial vertices.
            for (uint d : detector_array) {
                vertex_t v1 = std::make_pair(d, true);
                vertex_t v2 = std::make_pair(d, false);

                if (bp_matching.count(v1)) {
                    label_table[v1] = FREE;
                } else {
#ifdef BDC_DEBUG
                    std::cout << "\t" << v1.first << " is an S-Boy.\n";
#endif
                    label_table[v1] = SBOY;
                    vertices.push_back(v1);
                    visited.insert(v1);
                    prev[v1] = v1;
                }
                label_table[v2] = FREE;
            }
        } else {
            for (uint d : detector_array) {
                vertex_t v = std::make_pair(d, true);
                if (label_table[v] == SBOY) {
                    vertices.push_back(v);
                }
            }
        }
        new_stage = false;

        // Now, perform DFS and search for a free SINGLE girl.
        bool found = false;
        vertex_t end;
#ifdef BDC_DEBUG
        std::cout << "BFS:\n";
#endif
        while (!vertices.empty()) {
            vertex_t v1 = vertices.front();
#ifdef BDC_DEBUG
            std::cout << "\tPopped " << v1.first << "(" << v1.second << ")\n";
#endif
            vertices.pop_front();
            if (label_table[v1] == TGIRL && !bp_matching.count(v1)) {
                found = true;
                end = v1;
                break;
            }
            // Otherwise, traverse to neighbors.
            for (uint d : detector_array) {
                if (d == v1.first) {
                    continue;
                }
                vertex_t v2 = std::make_pair(d, !v1.second);
                if (visited.count(v2)) {
                    continue;
                }
                // Compute slack.
                std::pair<uint, uint> di_dj = std::make_pair(v1.first, v2.first);
                fp_t raw_weight = path_table[di_dj].distance;
                wgt_t w = (wgt_t) (raw_weight * MWPM_INTEGER_SCALE);
                wgt_t slack = w - u_table[v1] - u_table[v2];
                // We only go across edges with 0 slack.
                if (slack != 0) {
                    continue;
                }
#ifdef BDC_DEBUG
                std::cout << "\t\tWalking to " << v2.first 
                    << "(" << v2.second << ")\n";
#endif
                // Otherwise, try to label v2.
                // Rules:
                // (1) If v1 is SBOY and v2 is FREE and (v1, v2) is not
                //      in the matching, label v2 as TGIRL.
                // (2) If v2 is TGIRL and v1 is FREE and (v1, v2) is in
                //      the matching, label v1 as SBOY.
                if (label_table[v1] == SBOY 
                        && label_table[v2] == FREE
                        && (!bp_matching.count(v1) || bp_matching[v1] != v2))
                {
                    label_table[v2] = TGIRL;
#ifdef BDC_DEBUG
                    std::cout << "\t\tLabeled v2 as TGIRL.\n";
#endif 
                }
                if (label_table[v1] == TGIRL
                        && label_table[v2] == FREE
                        && bp_matching.count(v1)
                        && bp_matching[v1] == v2)
                {
                    label_table[v2] = SBOY;
#ifdef BDC_DEBUG
                    std::cout << "\t\tLabeled v2 as SBOY.\n";
#endif 
                }
                if (label_table[v2] == FREE) {
                    continue;
                }
                prev[v2] = v1;
                vertices.push_back(v2);
                visited.insert(v2);
            }
        }
        // If we succeed in finding a TGIRL, then use prev to 
        // augment the matching.
        if (found) {
#ifdef BDC_DEBUG
            std::cout << "Found augmenting path.\n";

            std::cout << "\t(before)Matching is size " 
                << bp_matching.size() << ".\n";
            for (auto mate : bp_matching) {
                vertex_t v1 = mate.first;
                vertex_t v2 = mate.second;
                std::cout << "\t\t" << v1.first << "("
                            << v1.second << ") --> " 
                            << v2.first <<  "(" 
                            << v2.second << ")\n";
            }
#endif
#ifdef BDC_DEBUG
            std::cout << "\tpath: " << end.first << "(" << end.second << ")";
#endif
            vertex_t curr = end;
            bool add = true;
            while (prev[curr] != curr) {
                vertex_t p = prev[curr];
#ifdef BDC_DEBUG
                std::cout << " " << p.first << "(" << p.second << ")";
#endif
                if (add) {
                    bp_matching[curr] = p;
                    bp_matching[p] = curr;
                }
                add = !add;
                curr = p;
            }
            new_stage = true;
#ifdef BDC_DEBUG
            std::cout << "\n";
            std::cout << "\t(after)Matching is size " 
                << bp_matching.size() << ".\n";
            for (auto mate : bp_matching) {
                vertex_t v1 = mate.first;
                vertex_t v2 = mate.second;
                std::cout << "\t\t" << v1.first << "("
                            << v1.second << ") --> " 
                            << v2.first <<  "(" 
                            << v2.second << ")\n";
            }
#endif
        } else {
            // Otherwise, update the slack variables.
            wgt_t slack_update = std::numeric_limits<wgt_t>::max();
#ifdef BDC_DEBUG
            std::cout << "Slack list:\n";
#endif
            for (uint di : detector_array) {
                vertex_t vi = std::make_pair(di, true);
                if (label_table[vi] != SBOY) {
                    continue;
                }
                if (u_table[vi] < slack_update) {
                    slack_update = u_table[vi];  
                }
                for (uint dj : detector_array) {
                    if (di == dj) {
                        continue;
                    }
                    vertex_t vj = std::make_pair(dj, false);
                    if (label_table[vj] != FREE) {
                        continue;
                    }

                    fp_t raw_weight = path_table[std::make_pair(di, dj)].distance;
                    wgt_t w = (wgt_t) (raw_weight * MWPM_INTEGER_SCALE);
                    wgt_t slack = w - u_table[vi] - u_table[vj];
#ifdef BDC_DEBUG
                    std::cout << "\tSlack of " << di << "," << dj << " = " 
                        << slack << "\n";
                    std::cout << "\t\tu1[" << di << "] = " << u_table[vi] << "\n";
                    std::cout << "\t\tu2[" << dj << "] = " << u_table[vj] << "\n";
#endif
                    if (slack < slack_update) {
                        slack_update = slack;
                    }
                }
            }

#ifdef BDC_DEBUG
            std::cout << "Made slack update: " << slack_update << "\n";
#endif
            for (uint d : detector_array) {
                vertex_t v1 = std::make_pair(d, true);
                if (label_table[v1] == SBOY) {
                    u_table[v1] += slack_update;
                }
                vertex_t v2 = std::make_pair(d, false);
                if (label_table[v2] == TGIRL) {
                    u_table[v2] -= slack_update;
                }
            }
            new_stage = true;
        }
    }
    // Finally, extract a matching from the bipartite matching.
    std::map<uint, uint> matching;
    for (auto mate : bp_matching) {
        vertex_t vi = mate.first;
        vertex_t vj = mate.second;
        if (vi.second) {
            matching[vi.first] = vj.first;
            matching[vj.first] = vi.first;
        }
    }
    return matching;
}

}   // gulliver
}   // qrc
