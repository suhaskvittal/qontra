/*
 *  author: Suhas Vittal
 *  date:   28 December 2023
 * */

#ifndef PROTEAN_VISUALIZATION_h
#define PROTEAN_VISUALIZATION_h

#include "qontra/protean/network.h"

#include <vtils/bijective_map.h>

#include <array>
#include <vector>
#include <string>

namespace qontra {
namespace protean {

// Useful typdefs:
template <size_t N> using coord_t=std::array<fp_t, N>;
typedef vtils::BijectiveMap<coord_t<2>, sptr<net::phys_vertex_t>> Plane;

// Main visualization functions:
//
// write_to_dot writes the PhysicalNetwork to an output_file.
struct render_config_t {
    Plane       plane;
    std::string layout_engine = "neato";

    bool do_not_render_out_of_plane_edges=false;
};

// Renders the entire network in one file.
void render_network(std::string output_file, PhysicalNetwork&, render_config_t);
/// Renders each check in a different file. Good for readability or fine-grained examination.
void render_network_by_check(std::string output_folder, std::string ext, PhysicalNetwork&, render_config_t);

// place_network maps each qubit in the PhysicalNetwork on a coordinate plane
// via linear programming. The program is a quadratic program and attempts
// to minimize the distance between qubits with links while preventing overlaps.
struct placement_config_t {
    fp_t    min_distance_between_qubits = 1.0;

    fp_t    edge_length_objective_coef = 1.0;
    fp_t    edge_crossing_objective_coef = 0.1;

    // Note that even though the graph is planar, a drawing may still have
    // edge crossings in practice (as it may be that a certain drawing requires
    // a really long edge).
    //
    // The edge crossing part of the LP formulation may introduce many
    // integer variables and constraints. Here are methods to reduce the impact
    // of the region:
    //
    // edge_crossing_skip: most drastic -- simply skip this part of the LP
    // edge_crossing_relax_variables: many of the slack variables are made continuous
    // edge_crossing_max_indicators: constraints and variables are only made until the
    //                              specified indicators have been made (or exit).
    bool    edge_crossing_skip = false;
    bool    edge_crossing_relax_variables = false;
    size_t  edge_crossing_max_indicators = 1'000'000;

    // Plane configuration:
    fp_t    x_min = 0.0;
    fp_t    x_max = 32.0;
    fp_t    y_min = 0.0;
    fp_t    y_max = 32.0;
};

Plane place_network(PhysicalNetwork&, placement_config_t);

}   // protean
}   // qontra

#endif  // PROTEAN_VISUALIZATION_h
