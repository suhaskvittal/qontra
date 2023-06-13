/*
 *  author: Suhas Vittal
 *  date:   7 June 2023
 * */

#ifndef DEPENDENCE_GRAPH_h
#define DEPENDENCE_GRAPH_h

#include "defs.h"
#Include "graph/graph.h"

namespace qontra {
namespace graph {

namespace dep {

template <class I_t>
struct vertex_t : base::vertex_t {
    I_t* inst_p;
};

struct edge_t : base::edge_t {
};

}   // dep

#define __DependenceGraphParent Graph<dep::vertex_t<I_t>, dep::edge_t>

// This is a class that represents a DAG for instructions.
// Two instructions have an edge in the DAG if they share an
// operand. The direction of the edge points in the direction
// of the later operation.
template <class I_t>
class DependenceGraph : public __DependenceGraphParent {
public:
    DependenceGraph(void)
        :__DependenceGraphParent(),
        root(new dep::vertex_t<I_t>),
        is_dag(true)
    {
        root->inst_p = nullptr;
    }

    DependenceGraph(const DependenceGraph& other)
        :__DependenceGraphParent(other),
        root(other.root),
        is_dag(other.is_dag)
    {}

    bool add_vertex(dep::vertex_t<I_t>* v) override {
        if (!Graph::add_vertex(v))  return false;

        std::map<uint, dag::vertex_t<I_t>*> operand_to_ancestor;
        for (auto x : v->inst_p->operands)  operand_to_ancestor[x] = root;
        // Get ancestors via BFS.
        search::callback_t<V_t> cb = [&] (V_t* v1, V_t* v2)
        {
            for (auto x : v2->inst_p->operands) operand_to_ancestor[x] = v2;
        };
        // Connect v to all ancestors
        std::set<dag::vertex_t<I_t>*> already_connected;
        for (auto x : v->inst_p->operands) {
            auto a = operand_to_ancestor[x];
            if (already_connected.count(a)) continue;

            dep::edge_t* e = new dep::edge_t;
            e->src = a;
            e->dst = v;
            e->is_undirected = false;
            Graph::add_edge(e);

            already_connected.insert(a);
        }
        return true;
    }

    bool add_edge(dep::edge_t* e) =delete;
    dag::vertex_t<I_t>* get_root(void) { return root; }
private:
    dep::vertex_t<I_t>* root;   // Single source of the DAG. Has no data (inst_p = nullptr).

    bool is_dag;
};

template <class I_t>
DependenceGraph<I_t> from_schedule(const schedule_t<I_t>& schedule) {
    DependenceGraph<I_t> graph;

    uint index = 0;
    for (auto& inst : schedule) {
        dep::vertex_t* v = new dep::vertex_t;
        v->id = index++;
        v->inst_p = &inst;
        graph.add_vertex(v);
    }
    return graph;
}

}   // graph
}   // qontra

#endif  // DEPENDENCE_GRAPH_h
