/* 
 *  author: Suhas Vittal
 *  date:   17 February 2024
 * */

#include "qontra/decoder/restriction.h"

#include <vtils/set_algebra.h>
#include <vtils/utility.h>

namespace qontra {

using namespace graph;
using namespace decoding;

void
push_back_assignment(std::vector<assign_t>& arr, const assign_t& m) {
    if (m.w == nullptr || m.v == m.w) {
        return;
    }
    arr.push_back(m);
}

void
erase_from_incidence_map(std::map<vpair_t, size_t>& incidence_map, const vpair_t& e) {
    if (incidence_map.count(e) && (--incidence_map[e]) == 0) {
        incidence_map.erase(e);
    }
}

void
update_correction(
        std::map<vpair_t, size_t>& incidence_map,
        stim::simd_bits_range_ref<SIMD_WIDTH> corr,
        fp_t& corr_log_pr,
        stim::simd_bits_range_ref<SIMD_WIDTH> local_corr,
        const std::set<vpair_t>& local_boundary,
        fp_t local_log_pr)
{
    corr ^= local_corr;
    for (const vpair_t& e : local_boundary) {
        erase_from_incidence_map(incidence_map, e);
    }
    corr_log_pr += local_log_pr;
}

bool
remove_widowed_edges(std::map<vpair_t, size_t>& incidence_map) {
    std::map<sptr<vertex_t>, size_t> vertex_inc_map;
    for (const auto& [e, cnt] : incidence_map) {
        const auto& [v1, v2] = e;
        vertex_inc_map[v1]++;
        vertex_inc_map[v2]++;
    }
    // Remove any pairs of vertices where both endpoints only have a single
    // incidence.
    bool any_removed = false;
    for (auto it = incidence_map.begin(); it != incidence_map.end(); ) {
        const auto& [v1, v2] = it->first;
        if (vertex_inc_map.at(v1) == 1 && vertex_inc_map.at(v2) == 1) {
            any_removed |= (it->second % 2 == 1);
            it = incidence_map.erase(it);
        } else {
            it++;
        }
    }
    return any_removed;
}

std::set<vpair_t>
find_face_subset_given_cc_map(
    const std::map<vpair_t, size_t>& x_cc_map,
    const std::set<face_t>& faces,
    stim::simd_bits_range_ref<SIMD_WIDTH> corr,
    fp_t& log_pr,
    std::vector<face_t>& applied_faces,
    sptr<vertex_t> vcen)
{
    log_pr = 0;
    // Go through faces and compute adjacency. We will do a BFS from a starting
    // face repeatedly to satisfy the boundary.
    struct v_t : base::vertex_t {
        const face_t* f_p;
    };
    Graph<v_t, base::edge_t> fgr;
    size_t index = 0;
    for (const face_t& fc : faces) {
        sptr<v_t> v = fgr.make_and_add_vertex(index++);
        v->f_p = &fc;
        // Make edges.
        for (sptr<v_t> w : fgr.get_vertices()) {
            if (v == w) continue;
            size_t cnt = 0;
            for (sptr<vertex_t> x : fc.vertices) {
                const auto& wfcv = w->f_p->vertices;
                if (std::find(wfcv.begin(), wfcv.end(), x) != wfcv.end()) {
                    cnt++;
                }
            }
            if (cnt == 2) fgr.make_and_add_edge(v,w);
        }
    }
    std::set<vpair_t> boundary_update;
    for (const auto& [vp, cnt] : x_cc_map) {
        if (boundary_update.count(vp)) continue;
        const auto& [x,y] = vp;
        if (x != vcen && y != vcen) continue;
        // Find two starting vertices with x and y in the face.
        std::vector<sptr<v_t>> starting;
        for (sptr<v_t> v : fgr.get_vertices()) {
            const auto& vfcv = v->f_p->vertices;
            if (std::find(vfcv.begin(), vfcv.end(), x) != vfcv.end()
                && std::find(vfcv.begin(), vfcv.end(), y) != vfcv.end())
            {
                starting.push_back(v);
            }
        }
        fp_t best_log_pr = std::numeric_limits<fp_t>::lowest();
        std::set<vpair_t> best_update;
        stim::simd_bits<SIMD_WIDTH> best_local_corr(corr.num_bits_padded());
        std::vector<face_t> best_applied_faces;
        for (sptr<v_t> s : starting) {
            // First check if face is wholly incident to x_cc_map.
            std::set<vpair_t> base_boundary;
            size_t s_inc_cnt = 0;
            for (size_t i = 0; i < s->f_p->vertices.size(); i++) {
                sptr<vertex_t> v = s->f_p->vertices[i];
                for (size_t j = i+1; j < s->f_p->vertices.size(); j++) {
                    sptr<vertex_t> w = s->f_p->vertices[j];
                    vpair_t p = make_vpair(v,w);
                    if (x_cc_map.count(p)) s_inc_cnt++;
                    base_boundary.insert(p);
                }
            }
            if (s_inc_cnt == 2) {
                best_log_pr = log(s->f_p->probability);
                best_update = base_boundary;
                best_local_corr.clear();
                best_applied_faces = { *(s->f_p) };
                for (uint64_t fr : s->f_p->frames) best_local_corr[fr] ^= 1;
                continue;
            }
            // We will move "left" and "right" from s, and choose the set of
            // faces that has the best probability.
            stim::simd_bits<SIMD_WIDTH> base_corr(corr.num_bits_padded());
            for (uint64_t fr : s->f_p->frames) base_corr[fr] ^= 1;
            for (sptr<v_t> ss : fgr.get_neighbors(s)) {
                fp_t local_log_pr = log(s->f_p->probability);
                std::set<vpair_t> local_boundary_update(base_boundary);
                stim::simd_bits<SIMD_WIDTH> local_corr(base_corr);
                std::vector<face_t> local_applied_faces{ *(s->f_p) };
                std::set<sptr<v_t>> visited{s};
                std::deque<sptr<v_t>> bfs{ss};
                while (bfs.size()) {
                    sptr<v_t> v = bfs.front();
                    bfs.pop_front();
                    if (visited.count(v)) continue;

                    for (uint64_t fr : v->f_p->frames) local_corr[fr] ^= 1;
                    local_log_pr += log(v->f_p->probability);
                    // Check if we have found something else in the boundary.
                    bool tried_to_remove_vp = false;
                    bool any_edge_in_cc_map = false;
                    std::set<vpair_t> tmp;
                    for (size_t i = 0; i < v->f_p->vertices.size(); i++) {
                        sptr<vertex_t> a = v->f_p->vertices[i];
                        for (size_t j = i+1; j < v->f_p->vertices.size(); j++) {
                            sptr<vertex_t> b = v->f_p->vertices[j];
//                          if (a != vcen && b != vcen) continue;
                            vpair_t p = make_vpair(a,b);
                            if (p == vp) {
                                // This is bad.
                                local_log_pr = std::numeric_limits<fp_t>::lowest();
                                tried_to_remove_vp = true;
                            }
                            if (x_cc_map.count(p) && !boundary_update.count(p)) {
                                any_edge_in_cc_map = true;
                            }
                            tmp.insert(p);
                        }
                    }
                    if (tried_to_remove_vp) break;
                    local_applied_faces.push_back(*(v->f_p));
                    local_boundary_update ^= tmp;
                    if (any_edge_in_cc_map) break;
                    for (sptr<v_t> w : fgr.get_neighbors(v)) {
                        bfs.push_back(w);
                    }
                    visited.insert(v);
                }
                if (local_log_pr > best_log_pr) {
                    best_log_pr = local_log_pr;
                    best_update = local_boundary_update;
                    best_local_corr = local_corr;
                    best_applied_faces = local_applied_faces;
                }
            }
        }
        boundary_update ^= best_update;
        corr ^= best_local_corr;
        log_pr += best_log_pr;
        vtils::push_back_range(applied_faces, best_applied_faces);
    }
    if (boundary_update.empty()) log_pr = std::numeric_limits<fp_t>::lowest();
    return boundary_update;
}

}   // qontra
