/*
 *  author: Suhas Vittal
 *  date:   28 December 2023
 * */

#ifndef PROTEAN_VISUALIZATION_h
#define PROTEAN_VISUALIZATION_h

#include "protean/network.h"

#include <vtils/bijective_map.h>

#include <array>
#include <vector>
#include <string>

namespace protean {

// Main visualization functions:
//
// write_to_dot writes the PhysicalNetwork to an output_file.
struct render_config_t {
    std::string layout_engine = "neato";

    bool do_not_render_out_of_plane_edges=false;
};

// Renders the entire network in one file.
void render_network(std::string output_file, PhysicalNetwork*, render_config_t);
/// Renders each check in a different file. Good for readability or fine-grained examination.
void render_network_by_check(std::string output_folder, std::string ext, PhysicalNetwork*, render_config_t);

// Attributes:
struct attr_list_t {
    std::string name;

    std::string style;

    std::string fillcolor;
    std::string fontcolor;
    std::string fontname;
    std::string fontsize;
    std::string shape;

    std::string headlabel;
    std::string taillabel;
    std::string penwidth;
    std::string weight;
};

attr_list_t get_attributes(sptr<net::phys_vertex_t>);
attr_list_t get_attributes(sptr<net::phys_edge_t>);

}   // protean

#endif  // PROTEAN_VISUALIZATION_h
