/*
 *  author: Suhas Vittal
 *  date:   31 December 2023
 * */

#ifndef GRAPH_ALGORITHMS_PLANARITY_h
#define GRAPH_ALGORITHMS_PLANARITY_h

#include "qontra/graph.h"

#include <limits>
#include <tuple>
#include <set>
#include <vector>

namespace qontra {
namespace graph {

template <class V, class E>
class LRPlanarity {
public:
    LRPlanarity(Graph<V, E>*);
    
    bool run(void);
private:
    // Useful typedefs:
    typedef std::tuple<sptr<V>, sptr<V>>    oriented_edge_t;

    struct interval_t {
        interval_t(sptr<E> x=nullptr, sptr<E> y=nullptr)
            :low(x),
            high(y)
        {}

        interval_t(const interval_t& other)
            :low(other.low),
            high(other.high)
        {}

        sptr<E> low;
        sptr<E> high;

        inline bool is_null() { return low == nullptr && high == nullptr; }

        inline bool operator==(const interval_t& other) const {
            return low == other.low && high == other.high;
        }

        inline bool operator!=(const interval_t& other) const {
            return !(*this == other);
        }
    };

    struct cpair_t {
        cpair_t(interval_t l=interval_t(), interval_t r=interval_t())
            :left(l),
            right(r)
        {}

        cpair_t(const interval_t& other)
            :left(other.left),
            right(other.right)
        {}

        interval_t left;
        interval_t right;

        inline bool is_null() { return left.is_null() && right.is_null(); }

        inline bool operator==(const cpair_t& other) const {
            return left == other.left && right == other.right;
        }

        inline bool operator!=(const cpair_t& other) const {
            return !(*this == other);
        }
    };

    void    r_dfs_1(sptr<V>);
    bool    r_dfs_2(sptr<V>);   // Returns false if nonplanarity found.

    // Helper functions:
    sptr<V> get_src(sptr<E>) const;
    sptr<V> get_dst(sptr<E>) const;

    bool    add_edge_constraints(sptr<E>, sptr<E> parent_edge);
    bool    is_conflicting(interval_t, sptr<E>);
    size_t  get_lowest(cpair_t);

    Graph<V, E>*            input_graph;
    std::vector<sptr<V>>    roots;

    std::map<sptr<V>, sptr<E>>          parent_edge_map;
    std::map<sptr<E>, oriented_edge_t>  orientation_map;

    std::map<sptr<V>, std::vector<sptr<E>>> outgoing_edge_map;
    // Phase-1 Variables:
    //
    // height: tree-path distance from root.
    std::map<sptr<V>, size_t>   height;
    std::map<sptr<E>, size_t>   lowpt;
    std::map<sptr<E>, size_t>   lowpt2;
    std::map<sptr<E>, size_t>   nesting_depth;
    // Phase-2 Variables:
    //
    //
    std::map<sptr<E>, sptr<E>>  ref;
    std::set<sptr<E>>           neg_side;
    std::vector<cpair_t>        conflict_pair_stack;
    std::map<sptr<E>, cpair_t>  stack_bottom;
    std::map<sptr<E>, sptr<E>>  lowpt_edge;

    template <class NUMBER> inline NUMBER inf() { return std::numeric_limits<NUMBER>::max(); }

    const interval_t NULL_INTERVAL;
    const cpair_t NULL_CONFLICT_PAIR;
};

template <class V, class E>
bool lr_planarity(Graph<V, E>*);

}   // graph
}   // qontra

#include "planarity.inl"

#endif  // GRAPH_ALGORITHMS_PLANARITY_h
