/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include "protean/visualization_attr.h"
#include "protean/visualization_lp.h"

#include <cxxviz.h>

#include <string.h>

namespace qontra {
namespace protean {

void
render_network(std::string output_file, PhysicalNetwork& network, render_config_t config) {
    GVC_t* gvc = gvContext();

    Agraph_t* gr = cxx_agopen("processor", Agundirected, NULL);

    std::map<sptr<net::phys_vertex_t>, Agnode_t*> phys_to_gvc;
    // Create a node for each vertex in network.
    for (sptr<net::phys_vertex_t> pv : network.get_vertices()) {
        attr_list_t attr = get_attributes(pv);
        Agnode_t* av = cxx_agnode(gr, attr.name, 1);
        phys_to_gvc[pv] = av;
        // Update the vertex attributes.
        cxx_agset(av, "fillcolor", attr.fillcolor);
        cxx_agset(av, "fontcolor", attr.fontcolor);
        cxx_agset(av, "fontname", attr.fontname);
        cxx_agset(av, "fontsize", attr.fontsize);
        cxx_agset(av, "label", attr.name);
        cxx_agset(av, "shape", attr.shape);
        cxx_agset(av, "style", attr.style);

        if (config.plane.count(pv)) {
            coord_t<2> p = config.plane.at(pv);
            std::string ps = std::to_string(p[0]) + "," + std::to_string(p[1]) + "!";
            cxx_agset(av, "pos", ps);
        }
    }
    // Create edges. Here, we will label the edges by the roles that have 
    // the edge.
    for (sptr<net::phys_edge_t> pe : network.get_edges()) {
        attr_list_t attr = get_attributes(pe);

        sptr<net::phys_vertex_t> pv = std::reinterpret_pointer_cast<net::phys_vertex_t>(pe->src),
                                 pw = std::reinterpret_pointer_cast<net::phys_vertex_t>(pe->dst);
        Agnode_t* av = phys_to_gvc[pv];
        Agnode_t* aw = phys_to_gvc[pw];
        Agedge_t* ae = cxx_agedge(gr, av, aw, attr.name, 1);

        cxx_agset(ae, "dir", "none");
        cxx_agset(ae, "labelfontname", attr.fontname);
        cxx_agset(ae, "labelfontsize", attr.fontsize);
        cxx_agset(ae, "headlabel", attr.headlabel);
        cxx_agset(ae, "penwidth", attr.penwidth);
        cxx_agset(ae, "style", attr.style);
        cxx_agset(ae, "taillabel", attr.taillabel);
        cxx_agset(ae, "weight", attr.weight);
    }
    // Get output_file extension.
    const char* filename = output_file.c_str();
    char* ext = strrchr(filename, '.') + 1; // +1 to omit the '.'
    // Render the graph.
    gvLayout(gvc, gr, config.layout_engine.c_str());
    gvRenderFilename(gvc, gr, ext, filename);
    gvFreeLayout(gvc, gr);
    agclose(gr);
}

Plane
place_network(PhysicalNetwork& network, placement_config_t config) {
    // Build the LP.
    LP mgr;
    lp_add_variables(mgr, network, config);
    lp_add_minimum_distance_constraints(mgr, network, config);

    lp_expr_t objective;
    if (!config.edge_crossing_skip) {
        objective += 
            config.edge_crossing_objective_coef * lp_add_crossing_edges_objective(mgr, network, config);
    }
    objective += config.edge_length_objective_coef * lp_add_edge_distance_objective(mgr, network);
    mgr.build(objective, false);
    // Now, we can get the result and build the plane.
    fp_t obj;
    int solstat;
    if (mgr.solve(&obj, &solstat)) {
        std::cerr << "place_network: program is infeasible" << std::endl;
        return Plane();
    }
    // Form layout from the LP results.
    Plane layout;
    for (sptr<net::phys_vertex_t> v : network.get_vertices()) {
        coord_t<2> p{ mgr.get_value(get_x(v)), mgr.get_value(get_y(v)) };
        layout.put(v, p);
    }
    return layout;
}

}   // protean
}   // qontra
