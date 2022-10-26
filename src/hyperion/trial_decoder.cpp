/*
 *  author: Suhas Vittal
 *  date:   2 October 2022
 * */

#include "hyperion/trial_decoder.h"
#include <limits>

namespace qrc {
namespace hyperion {

Blossom::Blossom(uint node_id)
    :node_id(node_id),
    is_blossom(false),
    label(Label::S),
    bfwdptrs(),
    bbwdptr(nullptr),
    is_visited(false),
    is_matched(false),
    mateptr(nullptr),
    prevptr(nullptr),
    ownerptr(this),
    u(0),
    z(0)
{}

bool
Blossom::contains(Blossom * b) {
    for (Blossom * x : bfwdptrs) {
        if (x == b) {
            return true;
        } else if (x->is_blossom && x->contains(b)) {
            return true; 
        }
    }
    return false;
}

void
Blossom::reset() {
    is_visited = false;
    prevptr = nullptr;
    ownerptr = this;
}

TrialDecoder::TrialDecoder(const stim::Circuit& circ)
    :MWPMDecoder(circ)
{}

DecoderShotResult
TrialDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
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
    auto matching = hungarian(detector_list);

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
TrialDecoder::hungarian(const std::vector<uint>& detector_array) {
    // We implement the Hungarian algorithm: https://en.wikipedia.org/wiki/Hungarian_algorithm
    typedef std::pair<uint, bool> vertex_t;

    std::map<vertex_t, vertex_t> bp_matching;
    std::map<vertex_t, wgt_t> potential;
    for (uint d : detector_array) {
        vertex_t v1 = std::make_pair(d, true);
        vertex_t v2 = std::make_pair(d, false);
        potential[v1] = 0;
        potential[v2] = 0;
    }

    while (bp_matching.size() < 2*detector_array.size()) {
        std::deque<vertex_t> vertex_queue;
        std::map<vertex_t, vertex_t> prev;
        std::set<vertex_t> visited;
        for (uint d : detector_array) {
            vertex_t v = std::make_pair(d, true);
            if (!bp_matching.count(v)) {
                vertex_queue.push_back(v);
                visited.insert(v);
                prev[v] = v;
            }
        }
        // Perform a BFS to find an augmenting path.
        bool found = false;
        vertex_t end;
        while (!vertex_queue.empty()) {
            vertex_t v1 = vertex_queue.back();
            vertex_queue.pop_back();
            // Found endpoint of augmenting path.
            if (!found && !bp_matching.count(v1) && !v1.second) {
                found = true;
                end = v1;
                continue;
            }
            for (uint d : detector_array) {
                if (d == v1.first) {
                    continue;
                }
                vertex_t v2 = std::make_pair(d, !v1.second);
                if (visited.count(v2)) {
                    continue;
                }
                if (v1.second && (bp_matching.count(v1) 
                    && bp_matching[v1] == v2)) 
                {
                    continue;
                }
                if (!v1.second && 
                    (!bp_matching.count(v1) || bp_matching[v1] != v2)) 
                {
                    continue;
                }
                // Check if edge is tight.
                auto di_dj = std::make_pair(v1.first, d);
                fp_t raw_weight = path_table[di_dj].distance;
                wgt_t w = (wgt_t)(raw_weight * MWPM_INTEGER_SCALE);
                if (potential[v1] + potential[v2] != w) {
                    continue;
                }
                prev[v2] = v1;
                visited.insert(v2);
                vertex_queue.push_back(v2);
            }
        }
        // If we found an augmenting path, then update
        // the matching. Otherwise, update the potentials.
        if (found) {
            vertex_t curr = end;
            bool add = true;
            while (prev[curr] != curr) {
                vertex_t p = prev[curr];
                if (add) {
                    bp_matching[p] = curr;
                    bp_matching[curr] = p;
                }
                add = !add;
                curr = p;
            }
        } else {
            // Compute potential update.
            wgt_t update = std::numeric_limits<wgt_t>::max();
            for (uint di : detector_array) {
                // We only consider edges where the source
                // was visited and the destination was not.
                vertex_t vi = std::make_pair(di, true);
                if (!visited.count(vi)) {
                    continue;
                }
                for (uint dj : detector_array) {
                    if (di == dj) {
                        continue;
                    }
                    vertex_t vj = std::make_pair(dj, false);
                    if (visited.count(vj)) {
                        continue;
                    }
                    auto di_dj = std::make_pair(di, dj);
                    fp_t raw_weight = path_table[di_dj].distance;
                    wgt_t w = (wgt_t)(raw_weight * MWPM_INTEGER_SCALE);
                    wgt_t slack = w - potential[vi] - potential[vj];
                    if (slack < update) {
                        update = slack;
                    }
                }
            }
            for (uint d : detector_array) {
                vertex_t v1 = std::make_pair(d, true);
                vertex_t v2 = std::make_pair(d, false);
                if (visited.count(v1)) {
                    potential[v1] += update;
                }
                if (visited.count(v2)) {
                    potential[v2] -= update;
                }
            }
        }
    }
    // Translate result into a matching of the original graph.
    std::map<uint, uint> matching;
    for (auto match : bp_matching) {
        vertex_t v1 = match.first;
        vertex_t v2 = match.second;
        matching[v1.first] = v2.first;
    }
    return matching;
}

std::map<uint, uint>
TrialDecoder::blossom(const std::vector<uint>& detector_array) {

    // Initialize a Blossom data structure for each vertex.
    std::map<uint, Blossom*> blossoms;
    for (uint d : detector_array) {
        Blossom * b = new Blossom(d);
        blossoms[d] = b;
    }
    // Start search.
    std::map<uint, uint> matching;
    while (matching.size() < detector_array.size()) {
        std::deque<Blossom*> queue;
        for (auto pair : blossoms) {
            Blossom * b = pair.second;
            if (b->bbwdptr) {   // This is a subblossom.
                continue;
            }
            if (b->is_matched) {
                b->label = Blossom::Label::Free;
            } else {
                b->label = Blossom::Label::S;
                queue.push_back(b);
            }
            b->reset();
        }

        Blossom * src1 = nullptr, *src2 = nullptr;
        bool found = false;
        while (!queue.empty()) {
            Blossom * b1 = queue.front();
            uint d1 = b1->node_id;
            queue.pop_front();
            for (auto pair : blossoms) {
                uint d2 = pair.first; 
                Blossom * b2 = pair.second;
                if (b2->is_visited) {
                    continue;
                }
                // Check if the edge has zero slack.
                // Compute weight.
                auto di_dj = std::make_pair(d1, d2);
                fp_t raw_weight = path_table[di_dj].distance;
                wgt_t w = (wgt_t)(raw_weight * MWPM_INTEGER_SCALE);
                // Compute z.
                wgt_t z = 0;
                for (auto pair : blossoms) {
                    Blossom * b = pair.second;
                    if (b->contains(b1) && b->contains(b2)) {
                        z += b->z;
                    }
                }
                wgt_t slack = w - b1->u - b2->u - z;
                if (slack != 0) {
                    continue;
                }


                if (b2->label == Blossom::Label::Free) {
                    // This has a mate. 
                    auto b3 = b2->mateptr;
                    b2->label = Blossom::Label::T;
                    b3->label = Blossom::Label::S;
                    b2->is_visited = true;
                    b3->is_visited = true;
                    queue.push_front(b3);
                    continue;
                } else if (b2->label == Blossom::Label::T) {
                    continue;
                }
                // If Label = S.
                if (b1->ownerptr != b2->ownerptr) {
                    src1 = b1;
                    src2 = b2;
                    found = true;
                    goto bfs_done;
                } 
                // Otherwise, we need to make a blossom.
                // First, find the root of blossom.
                Blossom * root;
                std::set<Blossom*> marked;
                Blossom * curr = b1;
                while (curr != nullptr) {
                    marked.insert(curr); 
                    curr = curr->prevptr;
                }
                curr = b2;
                while (curr != nullptr) {
                    if (marked.count(curr)) {
                        root = curr;
                        break;
                    }
                    curr = curr->prevptr;
                }
            }
        }
bfs_done:
        return matching;
    }
}

}   // hyperion
}   // qrc
