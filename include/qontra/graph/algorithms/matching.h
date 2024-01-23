/*
 *  author: Suhas Vittal
 *  date:   1 January 2024
 * */

#ifndef GRAPH_ALGORITHMS_MATCHING_h
#define GRAPH_ALGORITHMS_MATCHING_h

#include "qontra/graph.h"

#include <map>
#include <vector>

namespace qontra {
namespace graph {

// Matching Graph structures:
enum class label_t { s, t, free };

template <class V>
struct m_vertex_t : base::vertex_t {
    bool    is_blossom;
    label_t label;

    sptr<m_vertex_T<V>> label_src;

    // If is_blossom == true, then the following are well-defined:
    std::vector<sptr<V>>    cycle;
    sptr<V>                 base;
    // Otherwise, this is just a normal vertex:
    sptr<V> orig;
    sptr<V> containing_blossom;
};

struct m_edge_t : base::edge_t {};

template <class V> void rule_1(sptr<m_vertex_t<V>>, sptr<m_vertex_t<V>>);
template <class V> void rule_2(sptr<m_vertex_t<V>>, sptr<m_vertex_t<V>>);
template <class V> void rule_12(sptr<m_vertex_t<V>>, sptr<m_vertex_t<V>>);

template <class V> using augp_t=std::vector<sptr<m_vertex_t<V>>>;
template <class V> using cyc_t=std::vector<sptr<m_vertex_t<V>>>;

template <class T>
void cyclic_iterator_inc(T::iterator&, const T& container, bool forwards);

template <class V, class E>
class MatchingGraph : public Graph<m_vertex_t<V>, m_edge_t> {
public:
    MatchingGraph(Graph<V, E>*);

    // Backtracking result data structure:
    struct btr_result_t {
        augp_t<V>   augmenting_path_1;
        augp_t<V>   augmenting_path_2;

        cyc_t<V>            blossom_cycle;
        sptr<m_vertex_t<V>> blossom_base;

        std::vector<sptr<m_vertex_t<V>>> t_labeled_vertices;
    };

    void                assign_label(sptr<m_vertex_t<V>> to, sptr<m_vertex_t<V>> from, label_t label);

    btr_result_t        backtrack(sptr<m_vertex_t<V>>, sptr<m_vertex_t<V>>);
    void                apply_augmenting_path(augp_t<V>);
    sptr<m_vertex_t<V>> add_blossom(const std::vector<sptr<m_vertex_t<V>>>& cycle, sptr<m_vertex_t<V>> base);
    void                expand_blossom(sptr<m_vertex_t<V>>);

    void    rule_1(sptr<m_vertex_t<V>>, sptr<m_vertex_t<V>>);
    void    rule_2(sptr<m_vertex_t<V>>, sptr<m_vertex_t<V>>);
    void    rule_12(sptr<m_vertex_t<V>>, sptr<m_vertex_t<V>>);

    std::map<sptr<m_vertex_t<V>>, sptr<m_vertex_t<V>>> matching;
private:
    std::map<sptr<V>, sptr<m_vertex_t>> v_orig_to_m;

    const uint64_t BLOSSOM_ID_FLAG = 1l << 48;
};

template <class V, class E>
class MaxCardinalityMatching {
public:
    MaxCardinalityMatching(Graph<V, E>*);

    void stage(void);
private:
    std::deque<sptr<m_vertex_t<V>>> s_queue;

    Graph<V, E>*        input_graph;
    MatchingGraph<V, E> matching_graph;
};

template <class V, class E>
std::map<V, V> max_cardinality_matching(Graph<V, E>*);

}   // graph
}   // qontra

#include "matching.inl"

#endif  // GRAPH_ALGORITHMS_MATCHING_h
