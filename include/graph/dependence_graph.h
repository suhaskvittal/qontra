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

struct vertex_t : base::vertex_t {
    Instruction* inst_p;
};

struct edge_t : base::edge_t {
};

}   // dep

#define __DependenceGraphParent Graph<dep::vertex_t, dep::edge_t>

// This is a class that represents a DAG for instructions.
// Two instructions have an edge in the DAG if they share an
// operand. The direction of the edge points in the direction
// of the later operation.
class DependenceGraph : public __DependenceGraphParent {
public:
    DependenceGraph(void)
        :__DependenceGraphParent(),
        root(new dep::vertex_t),
        depth(0),
        vertex_to_depth(),
        barrier_index(0)
    {
        root->inst_p = nullptr;
        vertex_to_depth[root] = 0;
    }

    DependenceGraph(schedule_t& schedule) 
        :DependenceGraph()
    {
        uint index = 0;
        for (auto& inst : schedule) {
            auto v = new dep::vertex_t;
            v->id = index++;
            v->inst_p = &inst;
            add_vertex(v);
        }
    }

    DependenceGraph(const DependenceGraph& other)
        :__DependenceGraphParent(other),
        root(other.root),
        depth(other.depth),
        vertex_to_depth(other.vertex_to_depth),
        barrier_index(other.barrier_index)
    {}

    ~DependenceGraph(void) {
        if (this->dealloc_on_delete) {
            for (auto v : this->vertices) delete v->inst_p;
        }
    }

    bool add_vertex(dep::vertex_t* v) override {
        if (!__DependenceGraphParent::add_vertex(v))  return false;

        std::map<uint, dep::vertex_t*> operand_to_ancestor;
        for (uint x : v->inst_p->operands)  operand_to_ancestor[x] = root;
        // Get ancestors via BFS.
        search::callback_t<dep::vertex_t> cb = 
        [&] (dep::vertex_t* v1, dep::vertex_t* v2)
        {
            for (uint x : v2->inst_p->operands) operand_to_ancestor[x] = v2;
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

    void add_barrier(const std::vector<uint>& operands) {
        Instruction* barrier = new Instruction;
        barrier->name = "NOP";
        barrier->operands = operands; 
        auto vb = new dep::vertex_t;
        vb->id = (barrier_index++) | (1 << 31);
        vb->inst_p = barrier;
        add_vertex(vb);
    }

    // Adding edges manually is disabled to maintain the DAG property.
    bool add_edge(dep::edge_t*)     { return false; }
    void delete_edge(dep::edge_t*)  { }

    schedule_t to_schedule(void) {
        schedule_t sch;

        // Add the instructions through a BFS.
        std::set<dep::vertex_t*> visited;
        search::callback_t<dep::vertex_t> cb = 
        [&] (dep::vertex_t* v1, dep::vertex_t* v2)
        {
            if (visited.count(v2))  return;
            sch.push_back(*(v2->inst_p));
            visited.insert(v2);
        };
        search::xfs(this, root, cb, false);
        return sch;
    }

    std::vector<dep::vertex_t*> get_vertices_at_depth(uint d) {
        std::vector<dep::vertex_t*> layer;
        for (auto v : this->vertices) {
            if (vertex_to_depth[v] == d)    layer.push_back(v);
        }
        return layer;
    }

    dep::vertex_t*  get_root(void)      { return root; }
    uint            get_depth(void)     { return depth; }

    uint            get_depth_of(dep::vertex_t* v) { return vertex_to_depth[v]; }
private:
    dep::vertex_t* root;   // Single source of the DAG. Has no data (inst_p = nullptr).

    uint                            depth;
    std::map<dep::vertex_t*, uint>  vertex_to_depth;

    uint    barrier_index;
};

}   // graph
}   // qontra

#endif  // DEPENDENCE_GRAPH_h
