/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#include "qontra/graph/decoding_graph.h"

#include <math.h>

namespace qontra {
namespace graph {

std::vector<int>
get_complementary_colors_to(std::vector<int> clist, int number_of_colors) {
    std::set<int> clist_set(clist.begin(), clist.end());
    std::vector<int> compl_list;
    for (int c = 0; c < number_of_colors; c++) {
        if (!clist_set.count(c)) compl_list.push_back(c);
    }
    return compl_list;
}

uint64_t
get_color_boundary_index(int color) {
    return BOUNDARY_FLAG | (static_cast<uint64_t>(color+1) << 48);
}

fp_t
compute_weight(fp_t probability) {
    return -log(probability);
}

}   // graph
}   // qontra
