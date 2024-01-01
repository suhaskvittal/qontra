/*
 *  author: Suhas Vittal
 *  date:   1 January 2024
 * */

#include "graph_prebuilt.h"

#include <graph/graph.h>
#include <graph/algorithms/planarity.h>

#include <iostream>

#include <stdlib.h>

int main(int argc, char* argv[]) {
    int n = atoi(argv[1]);

    GRAPH kgr = make_complete_graph(n);
    if (lr_planarity(&kgr)) {
        return 0;
    }
    return 1;
}

