/*
 *  author: Suhas Vittal
 *  date:   28 December 2023
 * */

#include "protean/network.h"

#include <string>

namespace qontra {
namespace protean {

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
}   // qontra
