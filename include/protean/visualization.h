/*
 *  author: Suhas Vittal
 *  date:   28 December 2023
 * */

#ifndef PROTEAN_VISUALIZATION_h
#define PROTEAN_VISUALIZATION_h

#include "protean/network.h"

#include <defs/bijective_map.h>

#include <array>
#include <vector>
#include <string>

namespace qontra {
namespace protean {

// Useful typdefs:
template <size_t N> using coord_t=std::array<fp_t, N>;
typedef BijectiveMap<coord_t<2>, net::phys_vertex_t> Plane;

// Main visualization functions:
//
// write_to_dot writes the PhysicalNetwork to an output_file in the .dot format.
void write_to_dot(std::string output_file, PhysicalNetwork&);

// place_network maps each qubit in the PhysicalNetwork on a coordinate plane
// via linear programming. The program is a quadratic program and attempts
// to minimize the distance between qubits with links while preventing overlaps.
struct placement_config_t {
    double min_distance_between_qubits = 1.0;

    double x_min = 0.0;
    double x_max = 1000.0;
    double y_min = 0.0;
    double y_max = 1000.0;
};

Plane place_network(PhysicalNetwork&, placement_config_t);

}   // protean
}   // qontra

#endif  // PROTEAN_VISUALIZATION_h
