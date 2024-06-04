/*
 *  author: Suhas Vittal
 *  date:   31 May 2024
 * */

#include "placc/tree.h"

#include <deque>
#include <map>

namespace placc {
    
ShorTree
ShorTree::build_tree(
        FPN* fpn,
        sptr<fpn_v_t> head,
        const std::vector<sptr<fpn_v_t>>& leaves,
        const std::set<sptr<fpn_v_t>>& blocked_qubits)
{
    ShorTree tree;

    std::deque<sptr<fpn_v_t>> bfs{ head };
    std::map<sptr<fpn_v_t>, sptr<fpn_v_t>> prev_map;
    std::set<sptr<fpn_v_t>> visited;

    prev_map[head] = nullptr;
    while (bfs.size()) {
        sptr<fpn_v_t> v = bfs.front(); bfs.pop_front();
        if (visited.count(v)) continue;
        // Check if v is a leaf.
        if (std::find(leaves.begin(), leaves.end(), v) != leaves.end()) {
            // Then commit all qubits using the prev_map.
            sptr<fpn_v_t> prev = nullptr;
            sptr<fpn_v_t> curr = v;
            while (curr != nullptr) {
                if (!tree.contains(curr)) {
                    tree.add_vertex(curr);
                }
                if (prev != nullptr && !tree.contains(prev, curr)) {
                    tree.make_and_add_edge(prev, curr);
                }
                prev = curr;
                curr = prev_map.at(curr);
            }
        } else {
            for (sptr<fpn_v_t> w : fpn->get_neighbors(v)) {
                if (blocked_qubits.count(w)) continue;
                bfs.push_back(w);
                prev_map[w] = v;
            }
        }
        visited.insert(v);
    }
    tree.head = head;
    tree.leaves = leaves;
    return tree;
}

size_t
ShorTree::get_spacetime_cost() {
    // Depth computation: perform a BFS.
    size_t max_depth = 0;
    std::map<sptr<fpn_v_t>, size_t> depth_map;
    
    std::deque<sptr<fpn_v_t>> bfs{ head };
    std::set<sptr<fpn_v_t>> visited;
    depth_map[head] = 0;
    while (bfs.size()) {
        sptr<fpn_v_t> v = bfs.front(); bfs.pop_front();
        if (depth_map.at(v) > max_depth) {
            max_depth = depth_map.at(v);
        }
        for (sptr<fpn_v_t> w : get_neighbors(v)) {
            if (visited.count(w)) continue;
            depth_map[w] = depth_map.at(v)+1;
            bfs.push_back(w);
        }
        visited.insert(v);
    }
    // Width: count qubits with degree > 1.
    size_t width = 0;
    for (sptr<fpn_v_t> v : get_vertices()) {
        if (get_degree(v) > 1) width++;
    }
    return max_depth*width;
}

ShorTree::ShorTree()
    :Graph(),
    head(),
    leaves()
{}

}   // placc
