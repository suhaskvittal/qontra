/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include "protean/visualization.h"

#include <vtils/graphviz/cxxviz.h>

#include <string.h>

using namespace qontra;
using namespace graph;
using namespace vtils;

namespace protean {

using namespace net;

void
add_node_to_render_graph(
        Agraph_t* gr,
        sptr<phys_vertex_t> pv,
        render_config_t config) 
{
    attr_list_t attr = get_attributes(pv);
    Agnode_t* av = cxx_agnode(gr, "PQ" + print_v(pv), 1);
    // Update the vertex attributes.
    cxx_agset(av, "fillcolor", attr.fillcolor);
    cxx_agset(av, "fontcolor", attr.fontcolor);
    cxx_agset(av, "fontname", attr.fontname);
    cxx_agset(av, "fontsize", attr.fontsize);
    cxx_agset(av, "label", attr.name);
    cxx_agset(av, "shape", attr.shape);
    cxx_agset(av, "style", attr.style);
}

void
add_edge_to_render_graph(
        Agraph_t* gr, 
        sptr<phys_edge_t> pe,
        render_config_t config)
{
    if (config.do_not_render_out_of_plane_edges && pe->is_out_of_plane()) return;
    attr_list_t attr = get_attributes(pe);

    sptr<phys_vertex_t> pv = pe->get_source<phys_vertex_t>(),
                             pw = pe->get_target<phys_vertex_t>();
    Agnode_t* av = cxx_agnode(gr, "PQ" + print_v(pv), 0);
    Agnode_t* aw = cxx_agnode(gr, "PQ" + print_v(pw), 0);
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
render_network(std::string output_file, PhysicalNetwork* network, render_config_t config) {
    GVC_t* gvc = gvContext();

    Agraph_t* gr = cxx_agopen("processor", Agundirected, NULL);
    // Create a node for each vertex in network->
    for (sptr<phys_vertex_t> pv : network->get_vertices()) {
        add_node_to_render_graph(gr, pv, config);
    }
    // Create edges. Here, we will label the edges by the roles that have 
    // the edge.
    for (sptr<phys_edge_t> pe : network->get_edges()) {
        add_edge_to_render_graph(gr, pe, config);
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
        PhysicalNetwork* network,
        render_config_t config) 
{
    GVC_t* gvc = gvContext();

    uptr<RawNetwork>& raw_network = network->get_raw_connection_network();
    for (sptr<raw_vertex_t> rpq : raw_network->get_vertices()) {
        if (rpq->qubit_type != raw_vertex_t::type::xparity
                && rpq->qubit_type != raw_vertex_t::type::zparity)
        {
            continue;
        }
        auto& support = raw_network->get_support(rpq);
        // Draw all qubits in the support. In case any of them share physical qubits, track what
        // physical qubits have been drawn.
        Agraph_t* gr = cxx_agopen("processor", Agundirected, NULL);

        std::set<sptr<phys_vertex_t>> visited;
        for (sptr<raw_vertex_t> rv : support.all) {
            sptr<phys_vertex_t> pv = network->get_physical_qubit_for(rv);
            if (visited.count(pv)) continue;

            add_node_to_render_graph(gr, pv, config);
        }

        for (sptr<raw_vertex_t> rv : support.all) {
            sptr<phys_vertex_t> pv = network->get_physical_qubit_for(rv);
            for (sptr<raw_vertex_t> rw : support.all) {
                sptr<phys_vertex_t> pw = network->get_physical_qubit_for(rw);
                if (pv <= pw) continue;
                sptr<phys_edge_t> pe = network->get_edge(pv, pw);
                if (pe == nullptr) continue;

                add_edge_to_render_graph(gr, pe, config);
            }
        }
        // Write to file.
        std::string filename = output_folder + "/" + print_v(rpq) + "." + ext;
        const char* c_fname = filename.c_str();
        const char* c_ext = ext.c_str();

        gvLayout(gvc, gr, config.layout_engine.c_str());
        gvRenderFilename(gvc, gr, c_ext, c_fname);
        gvFreeLayout(gvc, gr);
        agclose(gr);
    }
    gvFreeContext(gvc);
}

attr_list_t
get_attributes(sptr<phys_vertex_t> pv) {
    attr_list_t attributes;

    attributes.name = print_v(pv);
    // Set remaining attributes.
    if (pv->has_role_of_type(raw_vertex_t::type::data)) {
        attributes.fillcolor = "#36454f";
        attributes.fontcolor = "white";
        attributes.shape = "circle";
        attributes.fontsize = "10";
    } else if (pv->has_role_of_type(raw_vertex_t::type::xparity)
                && pv->has_role_of_type(raw_vertex_t::type::zparity))
    {
        attributes.fillcolor = "#e65480";
        attributes.fontcolor = "black";
        attributes.shape = "square";
        attributes.fontsize = "14";
    } else if (pv->has_role_of_type(raw_vertex_t::type::xparity)) {
        attributes.fillcolor = "#bd2031";
        attributes.fontcolor = "white";
        attributes.shape = "square";
        attributes.fontsize = "14";
    } else if (pv->has_role_of_type(raw_vertex_t::type::zparity)) {
        attributes.fillcolor = "#4169e1";
        attributes.fontcolor = "black";
        attributes.shape = "square";
        attributes.fontsize = "14";
    } else if (pv->has_role_of_type(raw_vertex_t::type::flag)) {
        attributes.fillcolor = "#faf0e6";
        attributes.fontcolor = "black";
        attributes.shape = "octagon";
        attributes.fontsize = "10";
    } else {
        attributes.fillcolor = "#faf0e6";
        attributes.fontcolor = "red";
        attributes.shape = "triangle";
        attributes.fontsize = "10";
    }
    attributes.fontname = "serif";
    attributes.style = "filled";
    return attributes;
}

attr_list_t
get_attributes(sptr<phys_edge_t> pe) {
    attr_list_t attributes;
    // Get endpoints:
    sptr<phys_vertex_t> pv = std::reinterpret_pointer_cast<phys_vertex_t>(pe->src),
                        pw = std::reinterpret_pointer_cast<phys_vertex_t>(pe->dst);
    attributes.name = "e" + print_v(pv) + ":" + print_v(pw);
    // Now set the attributes.
    if (pe->is_out_of_plane()) {
        attributes.style = "dotted";
        attributes.fontname = "serif";
        attributes.fontsize = "8";
        attributes.headlabel = std::to_string(pe->tsv_layer);
        attributes.taillabel = attributes.headlabel;
        attributes.penwidth = "0.5";
        attributes.weight = "0";
    } else {
        attributes.style = "solid";
        attributes.penwidth = "2.0";
        attributes.weight = "10";
    }
    return attributes;
}

}   // protean
