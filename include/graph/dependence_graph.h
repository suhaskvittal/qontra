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
//
// For efficiency, we do not enforce the DAG property (only the
// directed property). However, we offer a method to check if
// the object is a DAG.
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

    // Enforce the directed property.
    bool add_edge(dep::edge_t* e) override { 
        if (e->is_undirected) return false;
        return __DependenceGraphParent::add_edge(e);
    }

    bool check_if_dag(void) { // Returns true if this graph is a DAG.
        update_state();
        return is_dag;
    }

    dag::vertex_t<I_t>* get_root(void) { return root; }
protected:
    bool update_state(void) override {
        if (!__DependenceGraphParent::update_state())   return false;        
        typedef dep::vertex_t<I_t> V_t;
        std::map<V_t, uint> level_map;
        search::callback_t<V_t> cb = [&] (V_t* v1, V_t* v2)
        {
            if (level_map.count(v2) && level_map[v2] < level_map[v1]) {
                this->is_dag = false;
            }
            level_map[v2] = level_map[v1]+1;
        };
        xfs(this, root, cb, false);
    }
private:
    dep::vertex_t<I_t>* root;   // Single source of the DAG. Has no data (inst_p = nullptr).

    bool is_dag;
};

template <class I_t>
DependenceGraph<I_t>    from_schedule(const schedule_t<I_t>&);

}   // graph
}   // qontra

#endif  // DEPENDENCE_GRAPH_h
