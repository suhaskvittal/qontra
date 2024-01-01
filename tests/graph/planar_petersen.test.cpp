/*
 *  author: Suhas Vittal
 *  date:   1 January 2024
 * */

#include "graph_prebuilt.h"

#include <graph/graph.h>
#include <graph/algorithms/planarity.h>

#include <stdlib.h>

int main() {
    GRAPH gr = fast_make_graph_with_k_vertices(10);
    for (size_t i = 0; i < 5; i++) {
        sptr<VERTEX> v1 = gr.get_vertex(i);
        // First make edge with i+1 mod 5 vertex.
        sptr<VERTEX> w1 = gr.get_vertex((i+1) % 5);
        gr.add_edge(fast_make_edge(v1, w1));
        // Now make inner edges.
        sptr<VERTEX> v2 = gr.get_vertex(i+5);
        gr.add_edge(fast_make_edge(v1, v2));
        sptr<VERTEX> w2 = gr.get_vertex(((i+2) % 5) + 5);
        gr.add_edge(fast_make_edge(v2, w2));
    }
    std::cout << "n = " << gr.n() << ", m = " << gr.m() << std::endl;
    return static_cast<int>(!lr_planarity(&gr));
}
