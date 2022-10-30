/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#include "mwpm_decoder.h"

namespace qrc {

MWPMDecoder::MWPMDecoder(const stim::Circuit& circ, uint max_detector) 
    :Decoder(circ),
    longest_error_chain(0),
    path_table(),
    max_detector(max_detector)
{
    path_table = compute_path_table(graph);
}

std::string
MWPMDecoder::name() {
    return std::string(MWPM_DECODER_NAME);
}

bool
MWPMDecoder::is_software() {
    return true;
}

uint64_t
MWPMDecoder::sram_cost() {
    // We examine the size of the path_table.
    uint64_t n_bytes_sram = 0;   
    for (auto kv_pair : path_table) {
        std::pair<uint, uint> di_dj = kv_pair.first;
        DijkstraResult res = kv_pair.second;
        // We assume that the di_dj pair is stored
        // in a hardware matrix. We don't count the 
        // SRAM required to store keys.
        // 
        // Instead, we will only consider the SRAM
        // required to store the entry.
        uint64_t bytes_for_path = sizeof(res.path[0]) * res.path.size();
        uint64_t bytes_for_distance = sizeof(res.distance);
        n_bytes_sram += bytes_for_path + bytes_for_distance;
    }
    return n_bytes_sram;
}

DecoderShotResult
MWPMDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
    // Build Boost graph for MWPM.
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();

    // Note to self: fault ids in pymatching are the frames in DecodingGraph.
    // Log start time.
#ifdef __APPLE__
    auto start_time = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
#else
    struct timespec start_time_data;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start_time_data);
    auto start_time = start_time_data.tv_nsec;
#endif
    // Count number of detectors.
    uint8_t syndrome_is_even = 0x1;
    std::vector<uint> detector_list;
    for (uint di = 0; di < n_detectors; di++) {
        auto syndrome_bit = syndrome[di];
        if (di > max_detector) {
            syndrome_bit = 0;
        }
        if (syndrome_bit) {
            syndrome_is_even ^= 0x1;
            detector_list.push_back(di);
        }
    }

    if (!syndrome_is_even) {
        // Add boundary to matching graph.
        detector_list.push_back(BOUNDARY_INDEX);
    }
    // Build Blossom V instance.
    uint n_vertices = detector_list.size();
    uint n_edges = n_vertices * (n_vertices + 1) / 2;  // Graph is complete.
    PerfectMatching pm(n_vertices, n_edges);
    pm.options.verbose = false;
    // Add edges.
    for (uint vi = 0; vi < n_vertices; vi++) {
        uint di = detector_list[vi];
        for (uint vj = vi + 1; vj < n_vertices; vj++) {
            uint dj = detector_list[vj];
            std::pair<uint, uint> di_dj = std::make_pair(di,dj);
            if (path_table[di_dj].distance >= 1000.0) 
            {
                continue;  // There is no path.
            }
            fp_t raw_weight = path_table[di_dj].distance;
            // Note that we have already typedef'd qfp_t as wgt_t.
            wgt_t edge_weight = (wgt_t) (MWPM_INTEGER_SCALE * raw_weight);
            pm.AddEdge(vi, vj, edge_weight);
        }
    }
    // Solve instance.
    pm.Solve();
    std::map<uint, uint> matching;
    for (uint vi = 0; vi < n_vertices; vi++) {
        uint vj = pm.GetMatch(vi);
        uint di = detector_list[vi];
        uint dj = detector_list[vj];
        // Update matching data structure.
        matching[di] = dj;
        matching[dj] = di;
    }
    // Compute logical correction.
    std::vector<uint8_t> correction = get_correction_from_matching(matching);
    // Stop time here.
#ifdef __APPLE__
    auto end_time = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
#else
    struct timespec end_time_data;
    clock_gettime(CLOCK_MONOTONIC_RAW, &end_time_data);
    auto end_time = end_time_data.tv_nsec;
#endif
    auto time_taken = end_time - start_time;
    // Build result.
    DecoderShotResult res = {
        time_taken,
        0.0, // TODO
        is_logical_error(correction, syndrome, n_detectors, n_observables),
        correction,
        matching
    };
    return res;
}

std::vector<uint8_t>
MWPMDecoder::get_correction_from_matching(const std::map<uint, uint>& matching) {
    std::set<uint> visited;
    std::vector<uint8_t> correction(circuit.count_observables(), 0);
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
            for (uint obs : edge.frames) {
                // Flip the bit.
                if (obs >= 0) {
                    correction[obs] = !correction[obs];
                }
            }
        }
        if (detector_path.size()-1 > longest_error_chain) {
            longest_error_chain = detector_path.size() - 1;
        }
        visited.insert(di);
        visited.insert(dj);
    }
    return correction;
}

}  // qrc
