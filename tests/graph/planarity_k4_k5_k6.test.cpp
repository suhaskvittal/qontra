/*
 *  author: Suhas Vittal
 *  date:   1 January 2024
 * */

#include "graph_prebuilt.h"

#include <graph/graph.h>
#include <graph/algorithms/planarity.h>

int main() {
    GRAPH k4 = make_complete_graph<4>(),
          k5 = make_complete_graph<5>(),
          k6 = make_complete_graph<6>();
    if (lr_planarity(k4) && !lr_planarity(k5) && !lr_planarity(k6)) {
        return 0;
    }
    return 1;
}

