/*
 *  author: Suhas Vittal
 *  date:   1 January 2024
 * */

#include <vtils/set_algebra.h>

#include <algorithm>
#include <deque>
#include <iterator>

namespace qontra {
namespace graph {

template <class T> inline void
cyclic_iterator_inc(T::iterator& it, const T& container, bool forwards) {
    if (forwards) {
        it++;
        if (it == container.end()) it = container.begin();
    } else {
        if (it == container.begin()) it = container.end();
        it--;
    }

}

template <class V, class E>
MatchingGraph<V, E>::MatchingGraph(Graph<V, E>* graph)
    :Graph(),
    matching(),
    v_orig_to_m()
{
    // Copy all vertices and edges from graph into here.
    for (sptr<V> v : graph->get_vertices()) {
        sptr<m_vertex_t<V>> mv = std::make_shared<m_vertex_t<V>>();
        mv->id = v->id;
        assign_label(mv, nullptr, label_t::free);
        mv->orig = v;

        add_vertex(mv);
        v_orig_to_m[v] = mv;
    }

    for (sptr<E> e : graph->get_edges()) {
        sptr<V> src = std::reinterpret_pointer_cast<V>(e->src),
                dst = std::reinterpret_pointer_cast<V>(e->dst);
        sptr<m_edge_t> me = std::make_shared<m_edge_t>();
        me->src = v_orig_to_m[src];
        me->dst = v_orig_to_m[dst];
        add_edge(me);
    }
}

template <class V, class E> MatchingGraph<V, E>::btr_result_t
MatchingGraph<V, E>::backtrack(sptr<m_vertex_t<V>> v, sptr<m_vertex_t<V>> w) {
    btr_result_t res;
    cyc_t<V> v_cycle, r_w_cycle;
    res.blossom_base = nullptr;

    sptr<m_vertex_t<V>> v_curr = v,
                        w_curr = w;
    
    const auto update = [&] (sptr<m_vertex_t<V>> curr, augp_t<V>& path, cyc_t<V>& cycle) 
    {
        if (curr == nullptr) return curr;

        path.push_back(curr);
        if (res.blossom_base == nullptr) {
            if (visited.count(curr)) {
                // Found a common vertex.
                res.blossom_base = curr;
            } else {
                visited += curr;
                cycle.push_back(curr);
            }
        }
        if (curr->label == label_t::t) {
            res.t_labeled_vertices.push_back(curr);
        }
        if (curr->containing_blossom == nullptr) return curr->label_src;
        // Now, if curr is part of a blossom, we need to be careful. We
        // want to take the even augmenting path to the root.
        sptr<m_vertex_t<V>> B = curr->containing_blossom;
        auto curr_where = std::find(B->cycle.begin(), B->cycle.end(), curr);
        auto base_where = std::find(B->cycle.begin(), B->cycle.end(), B->Base);
        if (curr_where == base_where) return curr->label_src;
        
        auto distance = std::distance(curr_where, base_where);
        // Now, we need to check if we go forward, or backward in the cycle.
        bool forwards = (distance > 0) ^ (distance % 2 == 0);
        // We just did curr_where, so increment it.
        cyclic_iterator_inc(curr_where, B->cycle, forwards);
        while (curr_where != base_where) {
            update(*curr_where, path, cycle);
            cyclic_iterator_inc(curr_where, B->cycle, forwards);
        }
        return B->base->label_src;
    };

    // Compute everything via controlled backtracking.
    std::set<sptr<m_vertex_t<V>>> visited;
    while (v_curr != nullptr || w_curr != nullptr) {
        v_curr = update(v_curr, res.augmenting_path_1, v_cycle);
        w_curr = update(w_curr, res.augmenting_path_2, r_w_cycle);
    }
    // Merge the two cycles:
    res.blossom_cycle.insert(res.blossom_cycle.end(), v_cycle.begin(), v_cycle.end());
    res.blossom_cycle.insert(res.blossom_cycle.end(), r_w_cycle.rbegin(), r_w_cycle.rend());
    return res;
}

template <class V, class E>
MatchingGraph<V, E>::apply_augmenting_path(augp_t<V> p) {
    for (size_t i = 1; i < p.size(); i++) {
        auto v = p[i-1],
             w = p[i];
        if (i % 2 == 1) {
            matching[v] = w;
            matching[w] = v;
        } // Edges where i % 2 == 0 will automatically be removed from matching.
    }
}

template <class V, class E> inline void
MatchingGraph<V, E>::rule_1(sptr<m_vertex_t<V>> v, sptr<m_vertex_t<V>> w) {
    if (matching[v] != w && v->label == label_t::s && w->label == label_t::free) {
        assign_label(w, v, label_t::t);
    }
}

template <class V, class E> inline void
MatchingGraph<V, E>::rule_2(sptr<m_vertex_t<V>> v, sptr<m_vertex_t<V>> w) {
    if (matching[v] == w && v->label == label_t::free && w->label == label_t::t) {
        assign_label(v, w, label_t::s);
    }
}

template <class V, class E> inline void
MatchingGraph<V, E>::rule_12(sptr<m_vertex_t<V>> v, sptr<m_vertex_t<V>> w) {
    if (matching[w] != w && v->label == label_t::s && w->label == label_t::free) {
        assign_label(w, v, label_t::t);
        sptr<m_vertex_t<V>> u = matching.at(w);
        assign_label(u, w, label_t::s);
    }
}

template <class V, class E>
MaxCardinalityMatching<V, E>::MaxCardinalityMatching(Graph<V, E>* graph)
    :s_queue(),
    input_graph(graph),
    matching_graph(graph)
{}

template <class V, class E> void
MaxCardinalityMatching<V, E>::stage() {
    // Traverse through the S-vertices.
    while (s_queue.size()) {
        sptr<m_vertex_t<V>> v = s_queue.front();
        s_queue.pop_front();

        if (v->is_shrunk) continue;

        if (v->label != label_t::s) continue;
        for (sptr<m_vertex_t<V>> w : matching_graph.get_neighbors()) {
            // Either w is free or itself is an S-vertex. We ignore when w is a T-vertex.
            if (w->label == label_t::t) continue;
            if (w->label == label_t::free) {
                // Apply R12 and move on.
                matching_graph.rule_12(v, w);
                continue;
            }
            // Otherwise, w is a S-vertex. Backtrack to the label sources.
            // If p1.back() != p2.back(), then all is good, just apply both augmenting paths.
            if (p1.back() != p2.back()) {
                // Note we need not merge the two paths, as the edge between v and w is an
                // even edge, and will automatically be removed by apply p1.
                matching_graph.apply_augmenting_path(p1);
                matching_graph.apply_augmenting_path(p2);
                continue;
            }
            // Otherwise, we've formed a cycle, and we need to make a blossom.
            std::vector<sptr<m_vertex_t<V>>> cycle(p1_cycle.begin(), p1_cycle.end());
            cycle.insert(cycle.end(), rp2_cycle.rbegin(), rp2_cycle.rend());

            auto blossom = matching_graph.add_blossom(cycle, base);
            matching_graph.assign_label(blossom, nullptr, label_t::s);
            
        }
    }
}

}   // graph
}   // qontra
