/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#include "decoding_graph.h"

namespace qrc {

DecodingGraph::DecodingGraph() {
    // Initialize data structures.
    base = decoding_graph_base();
    detector_to_index = std::map<uint, uint>();
    // Setup boundary node.
    boundary_coord = std::array<fp_t, N_COORD>();
    boundary_coord.fill((uint)-1);
    add_detector(BOUNDARY_INDEX, boundary_coord);
}

DecodingGraph::DecodingGraph(const DecodingGraph& other) {
    base = decoding_graph_base(other.base);
    detector_to_index = 
        std::map<uint, uint>(other.detector_to_index);
}

DecodingGraph::DecodingGraph(DecodingGraph&& other) {
    base = decoding_graph_base(other.base);
    detector_to_index = 
        std::map<uint, uint>(other.detector_to_index);
}

void
DecodingGraph::add_detector(uint det, std::array<fp_t, N_COORD>& coord) {
    if (detector_to_index.count(det)) {
        // Simple update the coord data.
        auto index = detector_to_index[det];
        base[index].coord = coord;
    } else {
        DecodingVertex vertex_data = {
            coord,
            det
        };
        auto index = boost::add_vertex(vertex_data, base);
        detector_to_index[det] = index;
    }
}

void
DecodingGraph::add_edge(uint det1, uint det2,
        fp_t weight, fp_t e_prob, std::set<uint>& frames) 
{
    auto v1 = detector_to_index[det1];
    auto v2 = detector_to_index[det2];
    
    DecodingEdge edge_data = {
        std::make_pair(det1, det2),
        weight,
        e_prob,
        frames
    };
    boost::add_edge(v1, v2, edge_data, base);
}

uint
DecodingGraph::get(uint det_id) {
    if (!detector_to_index.count(det_id)) {
        // Add it to the graph.
        add_detector(det_id, boundary_coord);
    }
    return detector_to_index[det_id];
}

DecodingGraph
to_decoding_graph(const stim::Circuit& qec_circ) {
    DecodingGraph graph;

    stim::DetectorErrorModel dem = 
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
            qec_circ,
            true,  // decompose_errors
            true,  // fold loops
            false, // allow gauge detectors
            1.0,   // approx disjoint errors threshold
            false, // ignore decomposition failures
            false
        );
    // Create callbacks.
    error_callback_f err_f = 
        [&graph](fp_t e_prob, std::vector<uint> dets,
                std::set<uint> frames) 
        {
            if (e_prob == 0 || dets.size() == 0 || dets.size() > 2) {
                return;  // Zero error probability -- not an edge.
            }
            
            if (dets.size() == 1) {
                // We are connecting to the boundary here.
                dets.push_back(BOUNDARY_INDEX);
            }
            // Now, we should only have two entries in det.
            uint det1 = dets[0];
            uint det2 = dets[1];
            // Get vertex indices in boost graph.
            uint vi1 = graph.get(det1);
            uint vi2 = graph.get(det2);
            // Returns a tuple (edge_index, does_exist)
            auto graph_edge = boost::edge(vi1, vi2, graph.base);
            if (graph_edge.second) {
                // Get old edge data.
                DecodingEdge old_edge_data = graph.base[graph_edge.first];
                fp_t old_e_prob = old_edge_data.error_probability;
                std::set<uint> old_frames = old_edge_data.frames;
                if (frames == old_frames) {
                    e_prob = 
                        e_prob * (1-old_e_prob) + old_e_prob * (1-e_prob);
                    // We will introduce a new edge index, so just
                    // delete this.
                    boost::remove_edge(graph_edge.first, graph.base);
                }
            }
            fp_t edge_weight = (fp_t)log((1-e_prob)/e_prob);
            graph.add_edge(det1, det2, edge_weight, e_prob, frames);
        };
    detector_callback_f det_f =
        [&graph](uint det, std::array<fp_t, N_COORD> coords) 
        {
            graph.add_detector(det, coords);
        };
    // Declare coord offset array.
    uint det_offset = 0;
    std::array<fp_t, N_COORD> coord_offset;
    coord_offset.fill(0);  // Zero initialize.
    // Use callbacks to build graph.
    _read_detector_error_model(dem, 1, det_offset, coord_offset,
                                err_f, det_f);
    return graph;
}

PathTable
compute_path_table(DecodingGraph& graph) {
    PathTable path_table;
    // Perform Dijkstra's algorithm on the graph.
    uint n_detectors = boost::num_vertices(graph.base);
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
            _update_path_table(path_table, graph, i, vj, distances, predecessors);
        }
        // Do the same for the boundary.
        _update_path_table(path_table, graph, i, BOUNDARY_INDEX, 
                distances, predecessors);
    }
    return path_table;
}

void 
_read_detector_error_model(
        const stim::DetectorErrorModel& dem, uint n_iter,
        uint& det_offset, std::array<fp_t, N_COORD>& coord_offset,
        error_callback_f err_f, detector_callback_f det_f) 
{
    while (n_iter--) {  // Need this to handle repeats.
        for (stim::DemInstruction inst : dem.instructions) {
            stim::DemInstructionType type = inst.type;
            if (type == stim::DemInstructionType::DEM_REPEAT_BLOCK) {
                // The targets for this instruction are
                // (1) number of repeats, and
                // (2) block number.
                uint n_repeats = (uint)inst.target_data[0].data;
                stim::DetectorErrorModel subblock = 
                    dem.blocks[inst.target_data[1].data];
                _read_detector_error_model(subblock, n_repeats,
                        det_offset, coord_offset, err_f, det_f);
            } else if (type == stim::DemInstructionType::DEM_ERROR) {
                std::vector<uint> detectors;
                std::set<uint> frames;
                 
                fp_t e_prob = (fp_t)inst.arg_data[0];
                for (stim::DemTarget target : inst.target_data) {
                    if (target.is_relative_detector_id()) {
                        // This is a detector, add it to the list.
                        detectors.push_back(
                                (uint)target.data + det_offset);
                    } else if (target.is_observable_id()) {
                        frames.insert(target.data); 
                    } else if (target.is_separator()) {
                        // This is just due to decomposition.
                        // Handle each part of the decomposition
                        // separately.
                        err_f(e_prob, detectors, frames);
                        // Clear detectors and frames.
                        // We have already done the callback.
                        detectors.clear();
                        frames.clear();
                    }
                }
                // Handle last error.
                err_f(e_prob, detectors, frames);
            } else if (type == stim::DemInstructionType::DEM_SHIFT_DETECTORS) {
                det_offset += inst.target_data[0].data;
                uint k = 0;
                for (double a : inst.arg_data) {
                    coord_offset[k++] += (fp_t)a;
                }
            } else if (type == stim::DemInstructionType::DEM_DETECTOR) {
                // Compute coordinates.
                std::array<fp_t, N_COORD> coords(coord_offset);
                uint k = 0;
                for (double a : inst.arg_data) {
                    coords[k++] += (fp_t)a; 
                }
                // Now go through all declared detectors.
                for (stim::DemTarget target : inst.target_data) {
                    det_f(target.data + det_offset, coords);
                }
            }
        }
    }
}

void
_update_path_table(PathTable& path_table,
        DecodingGraph& graph, uint src, uint dst,
        const std::vector<fp_t>& distances,
        const std::vector<uint>& predecessors)
{
    // Compute path.
    uint vi = graph.get(src);
    uint vj = graph.get(dst);
    uint curr = vj;
    fp_t distance = distances[vj];
    std::vector<uint> path;
    if (distance < 1000.0) {
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

};  // qrc
