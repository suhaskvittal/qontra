/*
 *  author: Suhas Vittal
 *  date:   31 December 2023
 * */

#include <vtils/set_algebra.h>

#include <algorithm>

namespace qontra {
namespace graph {

template <class V, class E>
LRPlanarity<V, E>::LRPlanarity(Graph<V, E>* graph)
    :input_graph(graph),
    roots(),
    parent_edge_map(),
    orientation_map(),
    outgoing_edge_map(),
    // Phase-1 Variables
    height(),
    lowpt(),
    lowpt2(),
    nesting_depth(),
    // Phase-2 Variables
    ref(),
    neg_side(),
    conflict_pair_stack(),
    stack_bottom(),
    lowpt_edge(),
    // Other variables.
    NULL_INTERVAL(nullptr, nullptr),
    NULL_CONFLICT_PAIR(NULL_INTERVAL, NULL_INTERVAL)
{
    // Initialize variables.
    for (sptr<V> v : input_graph->get_vertices()) {
        height[v] = inf<size_t>();
    }

    for (sptr<E> e : input_graph->get_edges()) {
        lowpt[e] = 0;
        lowpt2[e] = 0;
        nesting_depth[e] = 0;

        ref[e] = nullptr;
        stack_bottom[e] = cpair_t(NULL_INTERVAL, NULL_INTERVAL);
        lowpt_edge[e] = nullptr;
    }
}

template <class V, class E> bool
LRPlanarity<V, E>::run() {
    // Orient input graph
    for (sptr<V> v : input_graph->get_vertices()) {
        if (height[v] < inf<size_t>()) {
            continue;
        }
        height[v] = 0;
        roots.push_back(v);
        parent_edge_map[v] = nullptr;

        r_dfs_1(v);
    }

    for (sptr<V> r : roots) {
        if (!r_dfs_2(r)) return false; // non-planarity found.
    }
    return true;
}

template <class V, class E> void
LRPlanarity<V, E>::r_dfs_1(sptr<V> v) {
    sptr<E> parent_edge = parent_edge_map[v];
    for (sptr<V> w : input_graph->get_neighbors(v)) {
        sptr<E> e = input_graph->get_edge(v, w);
        if (orientation_map.count(e)) continue;
        // Orient e.
        orientation_map[e] = std::make_tuple(v, w);
        lowpt[e] = height[v];
        lowpt2[e] = height[v];

        outgoing_edge_map[v].push_back(e);

        if (height[w] == inf<size_t>()) {
            // This vertex has NOT been visited, implying e is a tree edge.
            parent_edge_map[w] = e;
            height[w] = height[v]+1;
            r_dfs_1(w);
        } else {
            // This is a back edge.
            lowpt[e] = height[w]; // w is closer to the root than v.
        }
        // Determine nesting depth.
        nesting_depth[e] = 2 * lowpt[e] + (lowpt2[e] < height[v]);
        // Update lowpoints of parent edge.
        if (parent_edge == nullptr) continue;
        if (lowpt[e] < lowpt[parent_edge]) {
            lowpt2[parent_edge] = std::min(lowpt[parent_edge], lowpt2[e]);
            lowpt[parent_edge] = lowpt[e];
        } else if (lowpt[e] > lowpt[parent_edge]) {
            lowpt2[parent_edge] = std::min(lowpt2[parent_edge], lowpt[e]);
        } else {
            lowpt2[parent_edge] = std::min(lowpt2[parent_edge], lowpt2[e]);
        }
    }
}

template <class V, class E> bool
LRPlanarity<V, E>::r_dfs_2(sptr<V> v) {
    sptr<E> parent_edge = parent_edge_map[v];

    if (outgoing_edge_map[v].empty()) return true;
    // Sort outgoing edges by nesting depth.
    std::sort(outgoing_edge_map[v].begin(), outgoing_edge_map[v].end(),
            [&] (auto e, auto f) {
                return nesting_depth.at(e) < nesting_depth.at(f);
            });

    sptr<E> first_outgoing_edge = outgoing_edge_map[v][0];
    for (sptr<E> e : outgoing_edge_map[v]) {
        if (conflict_pair_stack.empty()) {
            stack_bottom[e] = NULL_CONFLICT_PAIR;
        } else {
            stack_bottom[e] = conflict_pair_stack.back();
        }

        sptr<V> dst = get_dst(e);
        if (e == parent_edge_map[dst]) {
            // Then, this is a tree edge.
            if (!r_dfs_2(dst)) return false; // non-planarity detected.
        } else {
            lowpt_edge[e] = e;
            conflict_pair_stack.emplace_back(NULL_INTERVAL, interval_t(e, e));
        }
        // Intergrate new return edges.
        if (lowpt[e] < height[v]) {
            // Then, there is a return edge.
            if (e == first_outgoing_edge) {
                lowpt_edge[parent_edge] = lowpt_edge[first_outgoing_edge];
            } else {
                if (!add_edge_constraints(e, parent_edge)) {
                    return false; // non-planarity detected.
                }
            }
        }
    }
    // Remove back edges returning to parent.
    if (parent_edge == nullptr) return true; // v is a root.

    sptr<V> u = get_src(parent_edge);
    // Trim back edges ending at u.
    while (conflict_pair_stack.size() && get_lowest(conflict_pair_stack.back()) == height[u]) {
        cpair_t P = conflict_pair_stack.back();
        conflict_pair_stack.pop_back();
        if (P.left.low != nullptr) {
            neg_side += P.left.low;
        }
    }
    // Check if we have any more conflict pairs (at most one).
    if (conflict_pair_stack.size()) {
        cpair_t P = conflict_pair_stack.back();
        conflict_pair_stack.pop_back();
        // Trim left interval.
        while (P.left.high != nullptr && get_dst(P.left.high) == u) {
            P.left.high = ref[P.left.high];
        }
        if (P.left.high == nullptr && P.left.low != nullptr) {
            ref[P.left.low] = P.right.low;
            neg_side += P.left.low;
            P.left.low = nullptr;
        }
        // Trim right interval.
        while (P.right.high != nullptr && get_dst(P.right.high) == u) {
            P.right.high = ref[P.right.high];
        }
        if (P.right.high == nullptr && P.right.low != nullptr) {
            ref[P.right.low] = P.left.low;
            neg_side += P.right.low;
            P.right.low = nullptr;
        }
        conflict_pair_stack.push_back(P);
    }
    // Side of parent_edge is side of a highest return edge.
    if (lowpt[parent_edge] < height[u]) {
        // parent_edge has a return edge.
        cpair_t& cpx = conflict_pair_stack.back();
        sptr<E> h_left = cpx.left.high,
                h_right = cpx.right.high;
        if (h_left != nullptr && (h_right == nullptr || lowpt[h_left] > lowpt[h_right])) {
            ref[parent_edge] = h_left;
        } else {
            ref[parent_edge] = h_right;
        }
    }
    return true;
}

template <class V, class E> inline sptr<V> 
LRPlanarity<V, E>::get_src(sptr<E> e) const {
    auto oe = orientation_map.at(e);
    return std::get<0>(oe);
}

template <class V, class E> inline sptr<V>
LRPlanarity<V, E>::get_dst(sptr<E> e) const {
    auto oe = orientation_map.at(e);
    return std::get<1>(oe);
}

template <class V, class E> bool
LRPlanarity<V, E>::add_edge_constraints(sptr<E> e, sptr<E> parent_edge) {
    cpair_t P;
    // Merge return edges of e into P.R.
    while (conflict_pair_stack.back() != stack_bottom[e]) {
        cpair_t Q = conflict_pair_stack.back();
        conflict_pair_stack.pop_back();

        if (!Q.left.is_null()) std::swap(Q.left, Q.right);
        // Regardless the graph is nonplanar if both L and R are not null.
        if (!Q.left.is_null()) return false;
        // Merge intervals.
        if (lowpt[Q.right.low] > lowpt[parent_edge]) {
            if (P.right.is_null()) {
                P.right.high = Q.right.high;
            } else {
                ref[P.right.low] = Q.right.high;
            }
            P.right.low = Q.right.low;
        } else {
            ref[Q.right.low] = lowpt_edge[parent_edge];
        }
    }
    // Merge conflicting return edges into P.L
    while (is_conflicting(conflict_pair_stack.back().left, e)
            || is_conflicting(conflict_pair_stack.back().right, e))
    {
        cpair_t Q = conflict_pair_stack.back();
        conflict_pair_stack.pop_back();
        
        if (is_conflicting(Q.right, e)) std::swap(Q.left, Q.right);
        if (is_conflicting(Q.right, e)) return false;
        // Merge interval below lowpt[e] into P.R.
        ref[P.right.low] = Q.right.high;
        if (Q.right.low != nullptr) {
            P.right.low = Q.right.low;
        }
        if (P.left.is_null()) {
            // This is the topmost interval.
            P.left.high = Q.left.high;
        } else {
            ref[P.left.low] = Q.left.high;
        }
        P.left.low = Q.left.low;
    }
    // Push P if something was merged.
    if (!P.is_null()) conflict_pair_stack.push_back(P);
    return true;
}

template <class V, class E> inline bool
LRPlanarity<V, E>::is_conflicting(LRPlanarity<V, E>::interval_t i, sptr<E> e) {
    return !i.is_null() && lowpt[i.high] > lowpt[e];
}

template <class V, class E> inline size_t
LRPlanarity<V, E>::get_lowest(LRPlanarity<V, E>::cpair_t P) {
    if (P.left.is_null()) {
        return lowpt[P.right.low];
    } else if (P.right.is_null()) {
        return lowpt[P.left.low];
    } else {
        return std::min(lowpt[P.left.low], lowpt[P.right.low]);
    }
}

template <class V, class E> bool
lr_planarity(Graph<V, E>* graph) {
    // First, do simple Euler check.
    size_t n = graph->n(),
           m = graph->m();
    if (n > 2 && m > 3*n - 6) return false;

    LRPlanarity<V, E> lrp(graph);
    return lrp.run();
}

}   // graph
}   // qontra
