/*
 *  author: Suhas Vittal
 *  date:   31 May 2024
 * */

#ifndef PLACC_TREE_h
#define PLACC_TREE_h

#include "placc/fpn.h"

#include <deque>
#include <set>
#include <vector>

namespace qg=qontra::graph;

namespace placc {

class CXManager;

class ShorTree : public qg::Graph<fpn_v_t, fpn_e_t> {
public:
    ShorTree(ShorTree&&) =default;

    static ShorTree 
        build_tree(FPN*,
                    sptr<fpn_v_t> head,
                    const std::vector<sptr<fpn_v_t>>& leaves,
                    const std::set<sptr<fpn_v_t>>& blocked_qubits);

    void init_bfs(void);

    size_t get_spacetime_cost(void); // depth * width

    sptr<fpn_v_t> get_head(void) const;
    std::vector<sptr<fpn_v_t>> get_leaves(void) const;
private:
    ShorTree(void);

    sptr<fpn_v_t> head;
    std::vector<sptr<fpn_v_t>> leaves;
};

}   // placc

#include "inl/tree.inl"

#endif  // PLACC_TREE_h
