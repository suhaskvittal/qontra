/*
 *  author: Suhas Vittal
 *  date:   1 November 2022
 * */

#include "gulliver/trial_decoder.h"
#include <limits>

namespace qrc {
namespace gulliver {

TrialDecoder::TrialDecoder(const stim::Circuit& circuit,
        uint detectors_per_round, const std::filesystem::path& table_directory)
    :Decoder(circuit),
    init_error_table(),
    middle_error_table(),
    final_error_table(),
    detectors_per_round(detectors_per_round),
    max_candidates(1024),
    baseline(circuit),
    path_table(compute_path_table(graph))
{
    // Read files in directory.
    for (auto file : std::filesystem::directory_iterator(table_directory)) {
#ifdef GTR_DEBUG
        std::cout << "Reading from file: " << file.path() << "\n";
#endif
        std::ifstream in(file.path());
        load_from_file(in);
    }
}

std::string
TrialDecoder::name() {
    return "GulliverTRIAL";
}

bool
TrialDecoder::is_software() {
    return true;
}

DecoderShotResult
TrialDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
#ifdef GTR_DEBUG
    std::cout << "=======================================\n";
#endif
    struct CandidateCmp {
        bool operator()(const ErrorEvent& a, const ErrorEvent& b) {
            return a.second < b.second;
        }
    };

    const uint n_detectors = circuit.count_detectors();
    const uint n_observables = circuit.count_observables();
    const uint rounds = n_detectors/detectors_per_round;

    const uint mid_min_d = detectors_per_round + (detectors_per_round >> 1);
    const uint mid_max_d = mid_min_d + detectors_per_round;

    std::vector<ErrorEvent> candidates;
    std::map<ErrorEvent, std::set<DecodingGraph::Edge>> dependent_errors;
    std::map<ErrorEvent, bool> is_dependent;
    for (uint r = 0; r <= rounds; r++) {
        std::sort(candidates.begin(), candidates.end(), CandidateCmp());
        if (candidates.size() > max_candidates) {
            candidates = std::vector<ErrorEvent>(candidates.begin(),
                                                candidates.begin()+max_candidates);
        }

        uint min_d, max_d;
        if (r == 0) {
            min_d = 0;
            max_d = detectors_per_round >> 1;
        } else if (r == rounds) {
            min_d = n_detectors - (detectors_per_round >> 1);
            max_d = n_detectors;
        } else {
            min_d = r*detectors_per_round - (detectors_per_round >> 1);
            max_d = min_d + detectors_per_round;
        }

        std::vector<uint8_t> sub_syndrome(syndrome.begin() + min_d,
                                            syndrome.begin() + max_d);
#ifdef GTR_DEBUG
        std::cout << "[ Round " << r << " ]\n";
        std::cout << "\tsmin = " << min_d << ", smax = " << max_d << "\n";
        std::cout << "\tSyndrome:";
        for (uint8_t b : sub_syndrome) {
            std::cout << b+0;
        }
        std::cout << "\n";
#endif
        uint hw = std::accumulate(sub_syndrome.begin(),
                                    sub_syndrome.end(),
                                    0);
        if (hw == 0) {
            // Extend all dependent errors by a measurement error.
            std::vector<ErrorEvent> next_candidates;
            for (ErrorEvent event : candidates) {
                if (!is_dependent[event]) {
                    next_candidates.push_back(event);
                    continue;
                } else if (r == rounds) {
                    continue;
                }
                std::set<DecodingGraph::Edge> new_errors;
                std::set<DecodingGraph::Edge> new_dependents;
                fp_t weight = event.second;

                bool valid_extension = true;
                for (auto edge : event.first) {
                    new_errors.insert(edge);
                    if (dependent_errors[event].count(edge)) {
                        uint di = edge.detectors.first;
                        uint dj = edge.detectors.second;
                        uint dm;
                        if (dj == BOUNDARY_INDEX || di > dj) {
                            dm = di;
                        } else {
                            dm = dj;
                        }
                        uint ndm = graph.get_next_round(dm).detector;
                        auto extens = graph.get_edge(dm, ndm);
                        if (extens.id == -1) {
                            valid_extension = false;
                            break;
                        } else {
                            new_errors.insert(extens);
                            new_dependents.insert(extens);
                            weight += extens.edge_weight;
                        }
                    } 
                }
                // Only add the event if it has been extended appropriately.
                if (valid_extension) {
                    ErrorEvent new_event = std::make_pair(new_errors, weight);
                    next_candidates.push_back(new_event);
                    dependent_errors[new_event] = new_dependents;
                    is_dependent[new_event] = true;
                }
            }
            candidates = next_candidates;
            continue;
        }
        
        std::vector<ErrorEvent> events;
        if (r == 0) {
            if (!init_error_table.count(sub_syndrome)) {
                std::cout << "\tFailure!\n";
                break;  // We cannot progress.
            }
            events = init_error_table[sub_syndrome];
        } else if (r == rounds) {
            if (!final_error_table.count(sub_syndrome)) {
                std::cout << "\tFailure!\n";
                break;  // We cannot progress.
            }
            events = final_error_table[sub_syndrome];
        } else {
            if (!middle_error_table.count(sub_syndrome)) {
                std::cout << "\tFailure!\n";
                break;  // We cannot progress.
            }
            events = middle_error_table[sub_syndrome];
        }

        // Shift the events to the proper detectors.
        // Also filter out events if it is the first or last round.
        std::vector<ErrorEvent> shifted_events;
        if (r == 0 || r == rounds) {
            // No modification needed for the outer two detectors sets.
            shifted_events = events;
        } else {
            for (auto event : events) {
                std::set<DecodingGraph::Edge> shifted_edges;
                bool skip = false;
#ifdef GTR_DEBUG2
                std::cout << "\tEvent from table:\n";
#endif
                for (auto edge : event.first) {
                    uint di = edge.detectors.first;
                    uint dj = edge.detectors.second;
#ifdef GTR_DEBUG2
                    std::cout << "\t\t" << di << "_" << dj;
#endif
                    // Identify primary (in round) and secondary (out of round)
                    // detectors.
                    uint pd, sd;
                    if (di >= mid_min_d && di < mid_max_d) {
                        pd = di;
                        sd = dj;
                    } else {
                        pd = dj;
                        sd = di;
                    }
                    bool upwards = min_d >= mid_max_d;
                    while (pd < min_d || pd >= max_d) {
                        DecodingGraph::Vertex next_pv, next_sv;
                        if (upwards) {
                            next_pv = graph.get_next_round(pd);
                            next_sv = graph.get_next_round(sd);
                        } else {
                            next_pv = graph.get_prev_round(pd);
                            next_sv = graph.get_prev_round(sd);
                        }
                        if (next_pv.id == -1 || next_sv.id == -1) {
                            skip = true;
#ifdef GTR_DEBUG2
                            std::cout << " skipped: invalid id!\n";
#endif
                            goto shift_exit;
                        }
                        pd = next_pv.detector;
                        sd = next_sv.detector;
                    }
                    auto shifted_edge = graph.get_edge(pd, sd);
                    if (shifted_edge.id == -1) {
#ifdef GTR_DEBUG2
                        std::cout << " skipped: no edge exists!\n";
#endif
                        skip = true;
                        break;
                    } else {
#ifdef GTR_DEBUG2
                        std::cout << " shifted as " <<
                                    pd << "_" << sd << "\n";
#endif
                        shifted_edges.insert(shifted_edge);
                    }
                }
shift_exit:
                if (!shifted_edges.empty() && !skip) {
#ifdef GTR_DEBUG2
                    std::cout << "\tShifted event added.\n";
#endif
                    shifted_events.push_back(
                            std::make_pair(shifted_edges, event.second));
                }
            }
        }
        // Now, branch from existing candidates.
        std::vector<ErrorEvent> next_candidates;
#ifdef GTR_DEBUG
        std::cout << "\tParent events:\n";
        for (auto parent : candidates) {
            std::cout << "\t\t";
            if (is_dependent[parent]) {
                std::cout << "(D)";
            }
            std::cout << "w = " << parent.second << ":";
            for (auto edge : parent.first) {
                uint di = edge.detectors.first;
                uint dj = edge.detectors.second;
                std::cout << " " << di << "_" << dj;
            }
            std::cout << "\n";
        }
        std::cout << "\tTotal events: " << shifted_events.size() << "\n";
#endif
        for (auto event : shifted_events) {
            // Check dependencies:
            std::set<DecodingGraph::Edge> dependencies;
            // This event must intersect with an event in the 
            // previous round.
            bool e_is_depender = false;
            // Events in the next round must intersect with this
            // event.
            bool e_is_dependent = false;
            for (auto edge : event.first) {
                uint di = edge.detectors.first;
                uint dj = edge.detectors.second;
                if (di != BOUNDARY_INDEX) {
                    if (di >= max_d) {
                        e_is_dependent = true;
                        dependencies.insert(edge);
                    } else if (di < min_d) {
                        e_is_depender = true;
                    }
                }
                if (dj != BOUNDARY_INDEX) {
                    if (dj >= max_d) {
                        e_is_dependent = true;
                        dependencies.insert(edge);
                    } else if (dj < min_d) {
                        e_is_depender = true;
                    }
                }
            }
#ifdef GTR_DEBUG3
            std::cout << "\tEvent (depender = " << e_is_depender
                << ", dependent = " << e_is_dependent << "):";
            for (auto edge : event.first) {
                uint di = edge.detectors.first;
                uint dj = edge.detectors.second;
                std::cout << " " << di << "_" << dj;
                if (dependencies.count(edge)) {
                    std::cout << "(D)";
                }
            }
            std::cout << "\n";
#endif 
            if (candidates.size() == 0) {
                if (!e_is_depender) {
                    is_dependent[event] = e_is_dependent;
                    dependent_errors[event] = dependencies;
                    next_candidates.push_back(event);
#ifdef GTR_DEBUG3
                    std::cout << "\t\tadding event with weight " << event.second 
                        << "\n";
#endif
                }
                continue;
            }
            // Search for a parent to build from.
            fp_t min_weight = std::numeric_limits<fp_t>::max();
            ErrorEvent min_parent;
            for (auto parent : candidates) {
                fp_t w = event.second + parent.second;
                if (e_is_depender && is_dependent[parent]) {
                    // Check that all dependent errors are in
                    // the depender.
                    bool all_dependencies_found = true;
                    for (auto depending_edge : dependent_errors[parent]) {
                        if (!event.first.count(depending_edge)) {
                            all_dependencies_found = false;
                            break;
                        }
                        w -= depending_edge.edge_weight;
                    }
                    if (!all_dependencies_found) {
                        continue;
                    }
#ifdef GTR_DEBUG3
                    std::cout << "\t\tFound matching parent.\n";
#endif
                } else if (!e_is_depender && !is_dependent[parent]) {
                    // Continue onwards.
                } else {
                    // Do not add.
                    continue;
                }
                if (w < min_weight) {
                    min_weight = w;
                    min_parent = parent;
                }

#ifdef GTR_DEBUG3
                std::cout << "\t\tadding event with weight " << w << "\n";
#endif
                // Update event with min_parent.
                std::set<DecodingGraph::Edge> set_diff;
                std::set_symmetric_difference(event.first.begin(),
                                    event.first.end(),
                                    parent.first.begin(),
                                    parent.first.end(),
                                    std::inserter(set_diff, set_diff.begin()));
                ErrorEvent new_event;
                new_event.first = set_diff;
                new_event.second = w;
                is_dependent[new_event] = e_is_dependent;
                dependent_errors[new_event] = dependencies;
                next_candidates.push_back(new_event); 
            }
        }
        candidates = next_candidates;
#ifdef GTR_DEBUG
        std::cout << "\tCandidates: " << candidates.size() << "\n";
#endif
    }


    std::vector<uint8_t> correction(n_observables, 0);
    fp_t min_weight = std::numeric_limits<fp_t>::max();
    ErrorEvent best_candidate;
    for (auto candidate : candidates) {
        if (candidate.second < min_weight) {
            min_weight = candidate.second;
            best_candidate = candidate;
        }
    }
    // Get correction.
    for (auto edge : best_candidate.first) {
        auto true_edge = graph.get_edge(edge.detectors.first, 
                            edge.detectors.second);
#ifdef GTR_DEBUG
        std::cout << "edge: " << edge.detectors.first << ","
                    << edge.detectors.second << "\n";
        std::cout << "\tframes:\n";
#endif
        for (uint obs : true_edge.frames) {
            if (obs >= 0) {
                correction[obs] = !correction[obs];
            }
        } 
    }
#ifdef GTR_DEBUG
    if (candidates.empty()) {
        std::cout << "Empty error event list!\n";
    }
#endif

    bool is_error = is_logical_error(correction, syndrome, 
                        n_detectors, n_observables);
#ifdef GTR_DEBUG
    std::set<ErrorEvent> visited;
    std::cout << "Given correction: " << correction[0]+0 << "\n";
    std::cout << "Is logical error: " << is_error << "\n";
    if (is_error) {
        auto mwpm_res = baseline.decode_error(syndrome);
        std::cout << "MWPM had error = " << mwpm_res.is_logical_error << "\n";
        std::set<uint> visited;
        std::vector<uint8_t> mwpm_corr(circuit.count_observables(), 0);
        for (auto di_dj : mwpm_res.matching) {
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
                std::cout << "\tTraveling along " << wi << "_" << wj << "\n";
                auto edge = graph.get_edge(wi, wj);
                // The edge should exist.
                for (uint obs : edge.frames) {
                    // Flip the bit.
                    if (obs >= 0) {
                        mwpm_corr[obs] = !mwpm_corr[obs];
                    }
                }
            }
            visited.insert(di);
            visited.insert(dj);
        }
    }
    for (auto candidate : candidates) {
        if (visited.count(candidate)) {
            continue;
        }
        visited.insert(candidate);

        std::vector<uint8_t> gen_syndrome(n_detectors, 0);
        std::cout << "\tSpeculated errors (weight = " << candidate.second
            << "): {";
        for (auto edge : candidate.first) {
            std::cout << " " << edge.detectors.first << "_" 
                << edge.detectors.second; 
            if (edge.detectors.first != BOUNDARY_INDEX) {
                gen_syndrome[edge.detectors.first] ^= 1;
            }
            if (edge.detectors.second != BOUNDARY_INDEX) {
                gen_syndrome[edge.detectors.second] ^= 1;
            }
        }
        std::cout << " }\n";
        std::cout << "\tGenerated syndrome:";
        uint init_max_d = detectors_per_round >> 1;
        uint final_min_d = n_detectors - (detectors_per_round >> 1);
        std::cout << "\n\t";
        for (uint i = 0; i < init_max_d; i++) {
            std::cout << gen_syndrome[i]+0;
        }   
        for (uint i = init_max_d; i < final_min_d; i++) {
            if ((i - init_max_d) % detectors_per_round == 0) {
                std::cout << "\n\t";
            }
            std::cout << gen_syndrome[i]+0;
        }
        std::cout << "\n\t";
        for (uint i = final_min_d; i < n_detectors; i++) {
            std::cout << gen_syndrome[i]+0;
        }
        std::cout << "\n";
    }
#endif
    DecoderShotResult res = {
        0.0,
        0.0,
        is_error,
        correction,
        std::map<uint,uint>()
    };
    return res;
}

void
TrialDecoder::load_from_file(std::ifstream& in) {
    uint n_detectors = circuit.count_detectors();
    uint max_d_init = detectors_per_round >> 1;  // There are only one of 
                                                 // X or Z stabilizers in the
                                                 // first round.
    uint min_d_final = n_detectors - (detectors_per_round >> 1);   // The X (or Z)
                                                                   // stabilizers
                                                                   // are measured
                                                                   // to detect
                                                                   // measurement
                                                                   // errors.
    std::string syndrome_string;
    while (std::getline(in, syndrome_string, ',')) {
        std::string n_events_str;
        std::getline(in, n_events_str, ',');
        int n_events = std::stoi(n_events_str);

        std::set<DecodingGraph::Edge> error_set;
        bool is_init_error = false;
        bool is_final_error = false;
        while (n_events--) {
            std::string di_str, dj_str;
            std::getline(in, di_str, '_');
            std::getline(in, dj_str, ',');
            uint di = (uint)std::stoul(di_str),
                 dj = (uint)std::stoul(dj_str);
            auto edge = graph.get_edge(di, dj);
            if ((di < max_d_init && di != BOUNDARY_INDEX) 
                || (dj < max_d_init && dj != BOUNDARY_INDEX)) 
            {
                is_init_error = true;
            }
            if ((di >= min_d_final && di != BOUNDARY_INDEX)
                || (dj >= min_d_final && dj != BOUNDARY_INDEX))
            {
                is_final_error = true;
            }
            error_set.insert(edge);
        }

        std::string w_str;
        std::getline(in, w_str, '\n');
        fp_t w = std::stof(w_str);
        // Convert syndrome string to vector.
        std::vector<uint8_t> syndrome;
        for (uint i = 0; i < syndrome_string.size(); i++) {
            if (syndrome_string[i] == '1') {
                syndrome.push_back(1);
            } else {
                syndrome.push_back(0);
            }
        }

        if (is_init_error) {
            if (!init_error_table.count(syndrome)) {
                init_error_table[syndrome] = std::vector<ErrorEvent>();
            }
            init_error_table[syndrome].push_back(std::make_pair(error_set, w));
        } else if (is_final_error) {
            if (!final_error_table.count(syndrome)) {
                final_error_table[syndrome] = std::vector<ErrorEvent>();
            }
            final_error_table[syndrome].push_back(std::make_pair(error_set, w));
        } else {
            if (!middle_error_table.count(syndrome)) {
                middle_error_table[syndrome] = std::vector<ErrorEvent>();
            }
            middle_error_table[syndrome].push_back(std::make_pair(error_set, w));
        }
    }
}

}  // gulliver
}  // qrc
