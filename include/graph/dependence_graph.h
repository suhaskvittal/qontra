/*
 *  author: Suhas Vittal
 *  date:   7 June 2023
 * */

#ifndef DEPENDENCE_GRAPH_h
#define DEPENDENCE_GRAPH_h

#include "defs.h"
#include "graph/algorithms/search.h"
#include "graph/graph.h"
#include "instruction.h"

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
        depth(0),
        vertex_to_depth()
    {
        root->inst_p = nullptr;
        vertex_to_depth[root] = 0;
    }

    DependenceGraph(const schedule_t<I_t>& schedule) 
        :DependenceGraph()
    {
        uint index = 0;
        for (const auto& inst : schedule) {
            auto v = new dep::vertex_t<I_t>;
            v->id = index++;
            v->inst_p = &inst;
            add_vertex(v);
        }
    }

    DependenceGraph(const DependenceGraph& other)
        :__DependenceGraphParent(other),
        root(other.root),
        depth(other.depth),
        vertex_to_depth(other.vertex_to_depth)
    {}

    ~DependenceGraph(void) {
        if (this->dealloc_on_delete) {
            for (auto v : this->vertices) delete v->inst_p;
        }
    }

    bool add_vertex(dep::vertex_t<I_t>* v) override {
        if (!__DependenceGraphParent::add_vertex(v))  return false;

        std::map<uint, dep::vertex_t<I_t>*> operand_to_ancestor;
        for (uint x : v->inst_p->operands)  operand_to_ancestor[x] = root;
        // Get ancestors via BFS.
        search::callback_t<dep::vertex_t<I_t>> cb = 
        [&] (dep::vertex_t<I_t>* v1, dep::vertex_t<I_t>* v2)
        {
            if (vertex_to_depth[v1] > vertex_to_depth[v2])  std::cout << "ERROR ERROR ERROR\n";
            for (uint x : v2->inst_p->operands) operand_to_ancestor[x] = v2;
        };
        search::xfs(this, root, cb, false);
        // Connect v to all ancestors
        std::set<dep::vertex_t<I_t>*> already_connected;
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

    // Adding edges manually is disabled to maintain the DAG property.
    bool add_edge(dep::edge_t*)     { return false; }
    void delete_edge(dep::edge_t*)  { }

    schedule_t<I_t> to_schedule(void) {
        schedule_t<I_t> sch;

        // Add the instructions through a BFS.
        std::set<dep::vertex_t<I_t>*> visited;
        search::callback_t<dep::vertex_t<I_t>> cb = 
        [&] (dep::vertex_t<I_t>* v1, dep::vertex_t<I_t>* v2)
        {
            if (visited.count(v2))  return;
            sch.push_back(*(v2->inst_p));
            visited.insert(v2);
        };
        search::xfs(this, root, cb, false);
        return sch;
    }

    dep::vertex_t<I_t>* get_root(void)      { return root; }
    uint                get_depth(void)     { return depth; }

    uint                get_depth_of(dep::vertex_t<I_t>* v) { return vertex_to_depth[v]; }
private:
    dep::vertex_t<I_t>* root;   // Single source of the DAG. Has no data (inst_p = nullptr).

    uint                                depth;
    std::map<dep::vertex_t<I_t>*, uint> vertex_to_depth;
};

}   // graph
}   // qontra

#endif  // DEPENDENCE_GRAPH_h
