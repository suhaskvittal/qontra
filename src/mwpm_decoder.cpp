/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#include "mwpm_decoder.h"

MWPMDecoder::MWPMDecoder(const stim::Circuit& circ) 
:Decoder(circ), path_table()
{
    // Perform Dijkstra's algorithm on the graph.
    uint n_detectors = boost::num_vertices(graph.base);
    // Parallelize the Dijkstra's computation if the Openmp Flag
    // is on.
    for (uint i = 0; i < n_detectors; i++) {
        uint vi = graph.get(i);
        // Build data structures for call.
        auto index_map =
            boost::get(boost::vertex_index, graph.base);
        auto weight_map = 
            boost::get(&DecodingEdge::edge_weight, graph.base);
        std::vector<uint> predecessors(n_detectors);
        std::vector<fp_t> distances(n_detectors);

        boost::dijkstra_shortest_paths(
            graph.base,  // Input graph
            vi,     // Source vertex
            // Named parameters:
            predecessor_map(
                boost::make_iterator_property_map(
                    predecessors.begin(), index_map
                )
            ).distance_map(
                boost::make_iterator_property_map(
                    distances.begin(), index_map
                )
            ).weight_map(weight_map)
        );
        // Results should be in predecessors and distances.
        for (uint vj = i + 1; vj < n_detectors; vj++) {
            update_path_table(i, vj, distances, predecessors);
        }
        // Do the same for the boundary.
        update_path_table(i, BOUNDARY_INDEX, distances, predecessors);
    }
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
    auto start_time = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
    // Count number of detectors.
    uint8_t syndrome_is_even = 0x1;
    std::vector<uint> detector_list;
    for (uint di = 0; di < n_detectors; di++) {
        auto syndrome_bit = syndrome[di];
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
            if (path_table[di_dj].distance >= MAX_PRACTICAL_DISTANCE) 
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
    // Compute logical correction.
    std::vector<uint8_t> correction(n_observables, 0);
    std::set<uint> visited;
    std::map<uint, uint> matching;
    for (uint vi = 0; vi < n_vertices; vi++) {
        if (visited.count(vi)) {
            continue;
        }
        uint vj = pm.GetMatch(vi);
        uint di = detector_list[vi];
        if (di >= match_detectors_less_than) {
            continue;
        }
        uint dj = detector_list[vj];
        // Update matching data structure.
        matching[di] = dj;
        matching[dj] = di;
        std::pair<uint, uint> di_dj = std::make_pair(di, dj);
        // Check path between the two detectors.
        // This is examining the error chain.
        std::vector<uint> detector_path(path_table[di_dj].path);
        for (uint i = 1; i < detector_path.size(); i++) {
            // Get edge from decoding graph.
            auto wi = graph.get(detector_path[i-1]);
            auto wj = graph.get(detector_path[i]);
            auto edge = boost::edge(wi, wj, graph.base);
            // The edge should exist.
            for (uint obs : graph.base[edge.first].frames) {
                // Flip the bit.
                if (obs >= 0) {
                    correction[obs] = !correction[obs];
                }
            }
        }
        visited.insert(vi);
        visited.insert(vj);
    }
    // Stop time here.
    auto end_time = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
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

void
MWPMDecoder::update_path_table(uint src, uint dst,
        const std::vector<fp_t>& distances,
        const std::vector<uint>& predecessors)
{
    // Compute path.
    uint vi = graph.get(src);
    uint vj = graph.get(dst);
    uint curr = vj;
    fp_t distance = distances[vj];
    std::vector<uint> path;
    if (distance < MAX_PRACTICAL_DISTANCE) {
        while (curr != vi) {
            if (curr >= boost::num_vertices(graph.base) || path.size() > 100) {
                distance = std::numeric_limits<fp_t>::max();
                path.clear();
                goto failed_to_find_path;
            }
            uint det = graph.base[curr].detector;
            path.push_back(det);
            curr = predecessors[curr];
        }
        path.push_back(graph.base[curr].detector);
    }
failed_to_find_path:
    // Build result.
    DijkstraResult res = {path, distance};
    path_table[std::make_pair(src, dst)] = res;
    path_table[std::make_pair(dst, src)] = res;
}
