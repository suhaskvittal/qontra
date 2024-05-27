/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#include "qontra/graph/decoding_graph.h"
#include "qontra/graph/decoding_graph/edge_class.h"

#include <vtils/utility.h>

namespace qontra {
namespace graph {

using namespace decoding;

DecodingGraph::DecodingGraph(
        const DetailedStimCircuit& circuit, size_t flips_per_error, bool reweigh_for_detectors)
    :HyperGraph(),
    number_of_colors(circuit.number_of_colors_in_circuit),
    error_polynomial(),
    expected_errors(),
    dijkstra_graph_map(),
    flagged_dijkstra_graph_map(),
    distance_matrix_map(),
    flagged_distance_matrix_map(),
    base_probability_map(),
    active_flags(),
    flag_detectors(circuit.flag_detectors),
    edge_classes(),
    flag_class_map(),
    edge_class_map(),
    nod_edges(),
    all_edges(),
    flags_are_active(false),
    renorm_factor(1.0),
    reweigh_for_detectors(reweigh_for_detectors)
{
    stim::DetectorErrorModel dem =
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
                circuit,
                false,  // decompose errors
                true,   // fold loops
                false,  // allow gauge detectors (non-deterministic detectors)
                0.0,    // approx disjoint errors threshold
                false,  // ignore decomposition failures
                false
            );
    // Create boundaries first.
    if (number_of_colors == 0) {
        // Then we have a single boundary.
        uint64_t d = get_color_boundary_index(COLOR_ANY);
        sptr<vertex_t> vb = make_and_add_vertex(d);
        vb->is_boundary_vertex = true;
    } else {
        std::vector<sptr<vertex_t>> boundary_vertices;
        for (int c = 0; c < number_of_colors; c++) {
            uint64_t d = get_color_boundary_index(c);
            sptr<vertex_t> vb = make_and_add_vertex(d);
            vb->is_boundary_vertex = true;
            vb->color = c;
            boundary_vertices.push_back(vb);
        }
        // Make all of these boundaries connected to each other by a 1.0 probability hyperedge.
        sptr<hyperedge_t> e = make_and_add_edge(boundary_vertices);
        e->probability = 1.0;
    }
    // Do not add edges to the graph immediately. First, we will collect all of them, and then
    // analyze the edges and add them accordingly (see resolve_edges).
    std::vector<sptr<hyperedge_t>> tentative_edges;
    auto df = 
        [&] (uint64_t d)
        {
            if (!this->contains(d)) this->make_and_add_vertex(d);
        };
    auto ef =
        [&] (fp_t p, std::vector<uint64_t> detectors, std::set<uint64_t> frames)
        {
            if (p == 0) return;
            // Remove all flags from the detectors.
            std::vector<uint64_t> flags;
            for (auto it = detectors.begin(); it != detectors.end(); ) {
                if (flag_detectors.count(*it)) {
                    flags.push_back(*it);
                    it = detectors.erase(it);
                } else {
                    it++;
                }
            }
            if (detectors.size() == 0) {
                sptr<hyperedge_t> e = this->make_edge({});
                e->flags = std::set<uint64_t>(flags.begin(), flags.end());
                e->frames = frames;
                e->probability = p;
                nod_edges.push_back(e);
                all_edges.push_back(e);
                return;
            }

            std::vector<sptr<vertex_t>> vlist;
            for (uint64_t d : detectors) {
                if (!this->contains(d)) {
                    vlist.push_back(make_and_add_vertex_(d, circuit));
                } else {
                    vlist.push_back(this->get_vertex(d));
                }
            }
            sptr<hyperedge_t> e = make_edge(vlist);
            e->probability = p;
            e->frames = frames;
            e->flags = std::set<uint64_t>(flags.begin(), flags.end());
            tentative_edges.push_back(e);
        };
    size_t detector_offset = 0;
    read_detector_error_model(dem, 1, detector_offset, ef, df);

    resolve_edges(tentative_edges, flips_per_error);
    // Sort nod_edges by the number of flags (greatest to least). This is useful when
    // computing the renormalization factor.
    std::sort(nod_edges.begin(), nod_edges.end(),
            [&] (auto x, auto y)
            {
                return x->flags.size() > y->flags.size();
            });
}

inline sptr<vertex_t>
DecodingGraph::make_vertex(uint64_t id) const {
    sptr<vertex_t> v = HyperGraph::make_vertex(id);
    v->base = v;
    return v;
}

sptr<vertex_t>
DecodingGraph::make_and_add_vertex_(uint64_t d, const DetailedStimCircuit& circuit) {
    sptr<vertex_t> x = make_and_add_vertex(d);
    if (circuit.detector_color_map.count(d)) {
        x->color = circuit.detector_color_map.at(d);
    }
    if (circuit.detector_base_map.count(d)) {
        uint64_t bd = circuit.detector_base_map.at(d);
        if (d != bd && !this->contains(bd)) {
            sptr<vertex_t> _x = make_and_add_vertex_(bd, circuit);
            x->base = _x;
        } else {
            x->base = this->get_vertex(bd);
        }
    }
    return x;
}

void
DecodingGraph::resolve_edges(
        const std::vector<sptr<hyperedge_t>>& edge_list, size_t flips_per_error) 
{
    std::vector<EdgeClass> classes = EdgeClass::from_edges(edge_list);
    // Go through and add boundaries to any edges if necessary.
    for (EdgeClass& c : classes) {
        // Add boundaries based on whether or not the representative is a flag edge.
        sptr<hyperedge_t> rep = c.get_representative();
        if (rep->flags.empty()) {
            if (rep->get_order() < flips_per_error) { 
                if (number_of_colors == 0) {
                    if (rep->get_order() + 1 == flips_per_error) {
                        sptr<vertex_t> vb = get_boundary_vertex(COLOR_ANY);
                        c.add_vertex(vb);
                    }
                } else {
                    std::vector<sptr<vertex_t>> boundary_list = 
                        get_complementary_boundaries_to(rep->get<vertex_t>());
                    if (boundary_list.size() + rep->get_order() == flips_per_error) {
                        for (sptr<vertex_t> vb : boundary_list) c.add_vertex(vb);
                    }
                }
            }
        } else {
            // Here, we operate under the assumption that the representative edge is a two qubit error
            // that will flip two detectors.
            if (rep->get_order() == 1) {
                // We have a single detector that must be linked to a boundary. We can just link to the
                // boundary of the same color.
                int color = rep->get<vertex_t>(0)->color;
                c.add_vertex(get_boundary_vertex(color));
            }
        }

        // Update any edge probabilities.
        std::set<uint64_t> all_flags;
        for (sptr<hyperedge_t> e : c.get_edges()) {
            if (e->flags.empty()) {
                // Check if a similar edge already exists. If so, update that.
                sptr<hyperedge_t> _e = get_edge(e->get<vertex_t>());
                if (_e != nullptr) {
                    fp_t& p1 = e->probability,
                        & p2 = _e->probability;
                    if (e->frames == _e->frames) {
                        p2 = (1-p1)*p2 + (1-p2)*p1;
                    }
                } else {
                    add_edge(e);
                }
            }
            vtils::insert_range(all_flags, e->flags);
            edge_class_map[e] = c;
            all_edges.push_back(e);
        }
        edge_classes.push_back(c);
        for (uint64_t f : all_flags) {
            flag_class_map[f].push_back(c);
        }
    }
}

}   // graph
}   // qontra
