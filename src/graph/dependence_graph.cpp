/*
 *  author: Suhas Vittal
 *  date:   6 July 2023
 * */

#include "graph/dependence_graph.h"

namespace qontra {
namespace graph {

bool 
DependenceGraph::add_vertex(dep::vertex_t* v) {
    if (!__DependenceGraphParent::add_vertex(v))  return false;

    std::map<uint, dep::vertex_t*> operand_to_ancestor;
    for (uint x : v->inst_p->operands)  operand_to_ancestor[x] = root;
    // Get ancestors via BFS.
    search::callback_t<dep::vertex_t> cb = 
    [&] (dep::vertex_t* v1, dep::vertex_t* v2)
    {
        auto qubits = v2->inst_p->get_qubit_operands();
        for (uint x : qubits) operand_to_ancestor[x] = v2;
    };
    search::xfs(this, root, cb, false);
    // Connect v to all ancestors
    std::set<dep::vertex_t*> already_connected;
    uint max_ancestor_depth = 0;
    for (uint x : v->inst_p->operands) {
        auto a = operand_to_ancestor[x];
        if (already_connected.count(a)) continue;

        dep::edge_t* e = new dep::edge_t;
        e->src = a;
        e->dst = v;
        e->is_undirected = false;
        __DependenceGraphParent::add_edge(e);

        already_connected.insert(a);

        if (vertex_to_depth[a] > max_ancestor_depth) {
            max_ancestor_depth = vertex_to_depth[a];
        }
    }
    vertex_to_depth[v] = max_ancestor_depth+1;
    if (vertex_to_depth[v] > depth) depth = vertex_to_depth[v];
    return true;
}

}   // graph
}   // qontra
