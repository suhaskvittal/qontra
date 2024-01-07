/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

#include "qontra/protean/visualization_attr.h"

namespace qontra {
namespace protean {

using namespace net;

attr_list_t
get_attributes(sptr<phys_vertex_t> pv) {
    attr_list_t attributes;

    attributes.name = graph::print_v(pv);
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
        attributes.fillcolor = "black";
        attributes.fontcolor = "black";
        attributes.shape = "point";
        attributes.fontsize = "8";
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
    attributes.name = "e" + graph::print_v(pv) + ":" + graph::print_v(pw);
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
}   // qontra
