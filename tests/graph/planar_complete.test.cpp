/*
 *  author: Suhas Vittal
 *  date:   1 January 2024
 * */

#include "graph_prebuilt.h"

#include <graph/graph.h>
#include <graph/algorithms/planarity.h>

#include <stdlib.h>

int main(int argc, char* argv[]) {
    GRAPH gr;
    std::cout << "ARGC; " << argc << std::endl;
    if (argc == 2) {
        int n = atoi(argv[1]);
        gr = make_complete_graph(n);
    } else {
        int n1 = atoi(argv[1]),
            n2 = atoi(argv[2]);
        gr = make_bipartite_complete_graph(n1, n2);
    }
    return static_cast<int>(!lr_planarity(&gr));
}

