/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include "qontra/protean/visualization_attr.h"
#include "qontra/protean/visualization_lp.h"

#include <vtils/graphviz/cxxviz.h>

#include <string.h>

using namespace vtils;

namespace qontra {
namespace protean {

void
add_node_to_render_graph(
        Agraph_t* gr,
        PhysicalNetwork& net,
        sptr<net::phys_vertex_t> pv,
        render_config_t config) 
{
    attr_list_t attr = get_attributes(pv);
    Agnode_t* av = cxx_agnode(gr, "PQ" + graph::print_v(pv), 1);
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

        std::cout << "[ vis ] placing qubit " << graph::print_v(pv) << " at (" << p[0] << ", " << p[1] << ")\n";
        cxx_agset(av, "pos", ps);
    }
}

void
add_edge_to_render_graph(
        Agraph_t* gr, 
        PhysicalNetwork& network,
        sptr<net::phys_edge_t> pe,
        render_config_t config)
{
    if (config.do_not_render_out_of_plane_edges && pe->is_out_of_plane()) return;
    attr_list_t attr = get_attributes(pe);

    sptr<net::phys_vertex_t> pv = network.get_source(pe),
                             pw = network.get_target(pe);
    Agnode_t* av = cxx_agnode(gr, "PQ" + graph::print_v(pv), 0);
    Agnode_t* aw = cxx_agnode(gr, "PQ" + graph::print_v(pw), 0);
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

void
render_network(std::string output_file, PhysicalNetwork& network, render_config_t config) {
    GVC_t* gvc = gvContext();

    Agraph_t* gr = cxx_agopen("processor", Agundirected, NULL);
    // Create a node for each vertex in network.
    for (sptr<net::phys_vertex_t> pv : network.get_vertices()) {
        add_node_to_render_graph(gr, network, pv, config);
    }
    // Create edges. Here, we will label the edges by the roles that have 
    // the edge.
    for (sptr<net::phys_edge_t> pe : network.get_edges()) {
        add_edge_to_render_graph(gr, network, pe, config);
    }
    // Get output_file extension.
    const char* filename = output_file.c_str();
    char* ext = strrchr(filename, '.') + 1; // +1 to omit the '.'
    // Render the graph.
    gvLayout(gvc, gr, config.layout_engine.c_str());
    gvRenderFilename(gvc, gr, ext, filename);
    gvFreeLayout(gvc, gr);
    agclose(gr);
    gvFreeContext(gvc);
}

void
render_network_by_check(
        std::string output_folder,
        std::string ext,
        PhysicalNetwork& network,
        render_config_t config) 
{
    GVC_t* gvc = gvContext();

    RawNetwork raw_net = network.get_raw_connection_network();
    for (sptr<net::raw_vertex_t> rpq : raw_net.get_vertices()) {
        if (rpq->qubit_type != net::raw_vertex_t::type::xparity
                && rpq->qubit_type != net::raw_vertex_t::type::zparity)
        {
            continue;
        }
        auto& support = raw_net.get_support(rpq);
        // Draw all qubits in the support. In case any of them share physical qubits, track what
        // physical qubits have been drawn.
        Agraph_t* gr = cxx_agopen("processor", Agundirected, NULL);

        std::set<sptr<net::phys_vertex_t>> visited;
        for (sptr<net::raw_vertex_t> rv : support.all) {
            sptr<net::phys_vertex_t> pv = network.get_physical_qubit_for(rv);
            if (visited.count(pv)) continue;

            add_node_to_render_graph(gr, network, pv, config);
        }

        for (sptr<net::raw_vertex_t> rv : support.all) {
            sptr<net::phys_vertex_t> pv = network.get_physical_qubit_for(rv);
            for (sptr<net::raw_vertex_t> rw : support.all) {
                sptr<net::phys_vertex_t> pw = network.get_physical_qubit_for(rw);
                if (pv <= pw) continue;
                sptr<net::phys_edge_t> pe = network.get_edge(pv, pw);
                if (pe == nullptr) continue;

                add_edge_to_render_graph(gr, network, pe, config);
            }
        }
        // Write to file.
        std::string filename = output_folder + "/" + graph::print_v(rpq) + "." + ext;
        const char* c_fname = filename.c_str();
        const char* c_ext = ext.c_str();

        gvLayout(gvc, gr, config.layout_engine.c_str());
        gvRenderFilename(gvc, gr, c_ext, c_fname);
        gvFreeLayout(gvc, gr);
        agclose(gr);
    }
    gvFreeContext(gvc);
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
    int status;
    if ((status=mgr.solve(&obj, &solstat))) {
        std::cerr << "place_network: program is infeasible, solstat = " << solstat
                << ", status = " << status << std::endl;
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
