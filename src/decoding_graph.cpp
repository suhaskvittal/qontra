/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#include "decoding_graph.h"

namespace qrc {

static int32_t V_INDEX = 0;
static int32_t E_INDEX = 0;

static DecodingGraph::Vertex NULL_VERTEX(-1, std::array<fp_t, N_COORD>(), 0);
static DecodingGraph::Edge NULL_EDGE(-1, 0, 0, 0.0, 0.0, std::set<uint>());

DecodingGraph::DecodingGraph() 
    :detector_to_vertex(),
    vertices_to_edge(),
    vertex_to_next_round(),
    vertex_to_prev_round(),
    boundary_coord(),
    vertex_list(),
    adjacency_matrix()
{
    boundary_coord.fill((uint)-1);
    add_detector(BOUNDARY_INDEX, boundary_coord);
}

uint
DecodingGraph::count_detectors() {
    return vertex_list.size();
}

void
DecodingGraph::add_detector(uint det, std::array<fp_t, N_COORD>& coord) {
    if (detector_to_vertex.count(det) && detector_to_vertex[det] != NULL_VERTEX) {
        // Simple update the coord data.
        detector_to_vertex[det].coord = coord;
        for (auto& v : vertex_list) {
            if (v.detector == det) {
                v.coord = coord;
                break;
            }
        }
    } else {
        Vertex v(V_INDEX++, coord, det);
        detector_to_vertex[det] = v;
        vertex_list.push_back(v);
        adjacency_matrix[v] = std::vector<Vertex>();
    }
}

void
DecodingGraph::add_edge(uint det1, uint det2,
        fp_t weight, fp_t e_prob, std::set<uint>& frames) 
{
    Vertex& v1 = detector_to_vertex[det1];
    Vertex& v2 = detector_to_vertex[det2];
    
    Edge e(E_INDEX++, det1, det2, weight, e_prob, frames);
    vertices_to_edge[std::make_pair(v1,v2)] = e;
    vertices_to_edge[std::make_pair(v2,v1)] = e;
    adjacency_matrix[v1].push_back(v2);
    adjacency_matrix[v2].push_back(v1);
}

void
DecodingGraph::remove_vertex(const Vertex& v) {
    // Delete v from the vertex list
    for (auto it = vertex_list.begin(); it != vertex_list.end(); ) {
        if (*it == v) {
            it = vertex_list.erase(it); 
        } else {
            it++; 
        }
    }
    // And adjacency lists of its neighbors. Delete any
    // edges containing v.
    for (Vertex w : adjacency_matrix[v]) {
        auto adj_list = adjacency_matrix[w];
        for (auto it = adj_list.begin(); it != adj_list.end(); ) {
            if (*it == v) {
                it = adj_list.erase(it); 
                vertices_to_edge[std::make_pair(v,w)] = NULL_EDGE;
                vertices_to_edge[std::make_pair(w,v)] = NULL_EDGE;
            } else {
                it++; 
            }
        } 
    }
    adjacency_matrix[v].clear();
    detector_to_vertex[v.detector] = NULL_VERTEX;
}

void
DecodingGraph::remove_edge(const Edge& e) {
    Vertex v = get_vertex(e.detectors.first);
    Vertex w = get_vertex(e.detectors.second);
    vertices_to_edge[std::make_pair(v,w)] = NULL_EDGE;
    vertices_to_edge[std::make_pair(w,v)] = NULL_EDGE;
    // Remove v from adjacency list of w and vice versa.
    auto& adj_list_v = adjacency_matrix[v];
    auto& adj_list_w = adjacency_matrix[w];

    for (auto it = adj_list_v.begin(); it != adj_list_v.end(); ) {
        if (*it == w) {
            adj_list_v.erase(it); 
        } else {
            it++; 
        }
    }

    for (auto it = adj_list_w.begin(); it != adj_list_w.end(); ) {
        if (*it == v) {
            adj_list_w.erase(it); 
        } else {
            it++; 
        }
    }
}

DecodingGraph::Vertex
DecodingGraph::get_vertex(uint det_id) {
    if (!detector_to_vertex.count(det_id) 
            || detector_to_vertex[det_id] == NULL_VERTEX) 
    {
        // Add it to the graph.
        add_detector(det_id, boundary_coord);
    }
    return detector_to_vertex[det_id];
}

DecodingGraph::Vertex
DecodingGraph::get_next_round(uint det_id) {
    const Vertex v = get_vertex(det_id);
    return get_next_round(v);
}

DecodingGraph::Vertex
DecodingGraph::get_next_round(const Vertex& v) {
    if (vertex_to_next_round.count(v)) {
        return vertex_to_next_round[v];
    }
    if (v.detector == BOUNDARY_INDEX) {
        vertex_to_next_round[v] = v;
        return v;
    }

    Vertex next_round = NULL_VERTEX;
    for (Vertex w : vertex_list) {
        if (v.coord[0] == w.coord[0]
            && v.coord[1] == w.coord[1]
            && v.coord[2] + 1 == w.coord[2])
        {
            next_round = w;
            break;
        }
    }
    vertex_to_next_round[v] = next_round;
    return next_round;
}

DecodingGraph::Vertex
DecodingGraph::get_prev_round(uint det_id) {
    const Vertex& v = get_vertex(det_id);
    return get_prev_round(v);
}

DecodingGraph::Vertex
DecodingGraph::get_prev_round(const Vertex& v) {
    if (vertex_to_prev_round.count(v)) {
        return vertex_to_prev_round[v];
    }
    if (v.detector == BOUNDARY_INDEX) {
        vertex_to_prev_round[v] = v;
        return v;
    }

    Vertex prev_round = NULL_VERTEX;
    for (Vertex w : vertex_list) {
        if (v.coord[0] == w.coord[0]
            && v.coord[1] == w.coord[1]
            && v.coord[2] == w.coord[2] + 1)
        {
            prev_round = w;
            break;
        }
    }
    vertex_to_prev_round[v] = prev_round;
    return prev_round;
}

DecodingGraph::Edge
DecodingGraph::get_edge(uint det1, uint det2) {
    Vertex v1 = get_vertex(det1);
    Vertex v2 = get_vertex(det2);
    return get_edge(v1, v2);
}

DecodingGraph::Edge
DecodingGraph::get_edge(const Vertex& v1, const Vertex& v2) {
    auto v1_v2 = std::make_pair(v1, v2);
    auto v2_v1 = std::make_pair(v2, v1);
    if (vertices_to_edge.count(v1_v2)) {
        return vertices_to_edge[v1_v2]; 
    } else if (vertices_to_edge.count(v2_v1)) {
        return vertices_to_edge[v2_v1];
    } else {
        return NULL_EDGE; 
    }
}

uint32_t
DecodingGraph::get_chain_length(uint det1, uint det2) {
    Vertex src = get_vertex(det1);
    Vertex dst = get_vertex(det2);
    
    std::deque<Vertex> bfs_queue{src};
    std::set<Vertex> visited;
    std::map<Vertex, uint32_t> distance;
    distance[src] = 0;

    while (!bfs_queue.empty()) {
        Vertex v = bfs_queue.front();
        bfs_queue.pop_front();
        if (visited.count(v)) {
            continue;
        }

        for (Vertex w : adjacency_list(v)) {
//          if (w.detector != BOUNDARY_INDEX) {
                if (!distance.count(w) || distance[v] + 1 < distance[w]) {
                    distance[w] = distance[v] + 1;
                }
                bfs_queue.push_back(w);
//          }
        }
        visited.insert(v);
    }
    return distance[dst];
}

std::vector<DecodingGraph::Vertex>
DecodingGraph::vertices() {
    return vertex_list;
}

std::vector<DecodingGraph::Vertex>
DecodingGraph::adjacency_list(const Vertex& v) {
    return adjacency_matrix[v];
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
            auto graph_edge = graph.get_edge(det1, det2);
            if (graph_edge != NULL_EDGE) {
                // Get old edge data.
                fp_t old_e_prob = graph_edge.error_probability;
                std::set<uint> old_frames = graph_edge.frames;
                if (frames == old_frames) {
                    e_prob = 
                        e_prob * (1-old_e_prob) + old_e_prob * (1-e_prob);
                    // We will introduce a new edge index, so just
                    // delete this.
                    graph.remove_edge(graph_edge);
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
    auto vertices = graph.vertices();
    uint n_detectors = graph.count_detectors();
    for (uint i = 0; i < n_detectors; i++) {
        DecodingGraph::Vertex s = vertices[i];
        // Build data structures for call.
        std::map<DecodingGraph::Vertex, DecodingGraph::Vertex> predecessors;
        std::map<DecodingGraph::Vertex, fp_t> distances;

        _dijkstra(graph, s, distances, predecessors);

        for (uint j = i + 1; j < n_detectors; j++) {
            DecodingGraph::Vertex t = vertices[j];
            _update_path_table(path_table, graph, s, t, distances, predecessors);
        }
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
_dijkstra(DecodingGraph& graph, 
        const DecodingGraph::Vertex& src,
        std::map<DecodingGraph::Vertex, fp_t>& distances,
        std::map<DecodingGraph::Vertex, DecodingGraph::Vertex>& predecessors)
{
    typedef std::pair<DecodingGraph::Vertex, fp_t> PQVertex;

    struct DijkstraCmp {
        bool operator()(const PQVertex& v1, const PQVertex& v2) {
            return v1.second > v2.second; 
        };
    };
    std::map<DecodingGraph::Vertex, PQVertex> v2pv;
    std::priority_queue<PQVertex, std::vector<PQVertex>, DijkstraCmp> queue;
    for (DecodingGraph::Vertex v : graph.vertices()) {
        if (v == src) {
            distances[v] = 0; 
        } else {
            distances[v] = std::numeric_limits<fp_t>::max();
        } 
        predecessors[v] = v;

        PQVertex pv = std::make_pair(v, distances[v]);
        queue.push(pv);
        v2pv[v] = pv;
    }

    std::set<DecodingGraph::Vertex> visited;
    while (!queue.empty()) {
        PQVertex pvi = queue.top(); 
        auto vi = pvi.first;
        queue.pop();
        if (pvi.second != distances[vi]) {
            continue; 
        }

        auto adj_list = graph.adjacency_list(vi);
        for (DecodingGraph::Vertex vj : adj_list) {
            if (visited.count(vj)) {
                continue; 
            }
            DecodingGraph::Edge e = graph.get_edge(vi.detector, vj.detector);
            fp_t new_dist = distances[vi] + e.edge_weight;
            if (new_dist < distances[vj]) {
                distances[vj] = new_dist; 
                predecessors[vj] = vi;
                // Insert new entry into the priority queue.
                PQVertex pvj = std::make_pair(vj, new_dist);
                queue.push(pvj);
            }
        }
        visited.insert(vi);
    } 
}

void
_update_path_table(PathTable& path_table,
        DecodingGraph& graph, 
        const DecodingGraph::Vertex& src, 
        const DecodingGraph::Vertex& dst,
        std::map<DecodingGraph::Vertex, fp_t>& distances,
        std::map<DecodingGraph::Vertex, DecodingGraph::Vertex>& predecessors)
{
    // Compute path.
    DecodingGraph::Vertex curr = dst;
    fp_t distance = distances[dst];
    std::vector<uint> path;

    while (curr != src) {
        if (curr == predecessors[curr]) {
            distance = std::numeric_limits<fp_t>::max();
            path.clear();
            goto failed_to_find_path;
        }
        uint det = curr.detector;
        path.push_back(det);
        curr = predecessors[curr];
    }
    path.push_back(curr.detector);
failed_to_find_path:
    // Build result.
    DijkstraResult res = {path, distance};
    path_table[std::make_pair(src.detector, dst.detector)] = res;
    path_table[std::make_pair(dst.detector, src.detector)] = res;
}

}  // qrc
