/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

#include <math.h>

namespace qontra {
namespace graph {

inline sptr<decoding::vertex_t>
DecodingGraph::get_boundary_vertex(int color) {
    uint64_t db = get_color_boundary_index(color);
    return get_vertex(db);
}

inline error_chain_t
DecodingGraph::get(uint64_t d1, uint64_t d2) {
    return get(get_vertex(d1), get_vertex(d2));
}

inline error_chain_t
DecodingGraph::get(sptr<decoding::vertex_t> v1, sptr<decoding::vertex_t> v2) {
    return get(COLOR_ANY, COLOR_ANY, v1, v2);
}

inline error_chain_t
DecodingGraph::get(int c1, int c2, uint64_t d1, uint64_t d2) {
    return get(c1, c2, get_vertex(d1), get_vertex(d2));
}

inline error_chain_t
DecodingGraph::get(int c1, int c2, sptr<decoding::vertex_t> v1, sptr<decoding::vertex_t> v2) {
    update_state();
    if (c1 > c2) std::swap(c1, c2);
    auto c1_c2 = std::make_pair(c1, c2);
    auto& dm_map = flags_are_active ? flagged_distance_matrix_map[c1_c2] : distance_matrix_map[c1_c2];
    if (!dm_map.count(v1)) {
        dijkstra_(c1, c2, v1);
    }
    return dm_map[v1][v2];
}

template <class CONTAINER> std::vector<sptr<decoding::vertex_t>>
DecodingGraph::get_complementary_boundaries_to(CONTAINER vlist) {
    std::set<int> colors_in_vlist(vlist.begin(), vlist.end());
    for (sptr<decoding::vertex_t> v : vlist) colors_in_vlist.insert(v->color);

    std::vector<sptr<decoding::vertex_t>> boundary_list;
    for (int c : get_complementary_colors_to(colors_in_vlist, number_of_colors)) {
        boundary_list.push_back(get_boundary_vertex(c));
    }
    return boundary_list;
}

inline void
DecodingGraph::activate_flags(const std::vector<uint64_t>& all_detectors) {
    deactivate_flags();
    for (uint64_t d : all_detectors) {
        if (flag_detectors.count(d)) active_flags.push_back(d);
    }
    flags_are_active = !active_flags.empty();
}

inline void
DecodingGraph::deactivate_flags() {
    flagged_distance_matrix_map.clear();
    flagged_dijkstra_graph_map.clear();
    active_flags.clear();
    flags_are_active = false;
}

inline poly_t
DecodingGraph::get_error_polynomial() {
    update_state();
    return error_polynomial;
}

inline fp_t
DecodingGraph::get_expected_errors() {
    update_state();
    return expected_errors;
}

inline bool
DecodingGraph::update_state() {
    if (!HyperGraph::update_state()) return false;
    dijkstra_graph_map.clear();
    distance_matrix_map.clear();
    build_error_polynomial();
    return true;
}

template <class CONTAINER> inline std::vector<int>
get_complementary_color_to(CONTAINER clist, int number_of_colors) {
    std::set<int> clist_set(clist.begin(), clist.end());
    std::vector<int> compl_list;
    for (int c = 0; c < number_of_colors; c++) {
        if (!clist.count(c)) compl_list.push_back(c);
    }
    return compl_list;
}

inline uint64_t
get_color_boundary_index(int color) {
    return BOUNDARY_FLAG | (static_cast<uint64_t>(color+1) << 32);
}

inline fp_t
compute_weight(fp_t probability) {
    return -log(probability);
}

template <class F1, class F2> void
read_detector_error_model(
        const stim::DetectorErrorModel& dem,
        size_t n_iter,
        size_t& detector_offset,
        F1 ef,
        F2 df)
{
    while (n_iter--) {
        for (stim::DemInstruction inst : dem.instructions) {
            stim::DemInstructionType type = inst.type;
            if (type == stim::DemInstructionType::DEM_REPEAT_BLOCK) {
                size_t n_rep = static_cast<size_t>(inst.target_data[0].data);
                stim::DetectorErrorModel blk = dem.blocks[inst.target_data[1].data];
                read_detector_error_model(blk, n_rep, detector_offset, ef, df);
            } else if (type == stim::DemInstructionType::DEM_ERROR) {
                std::vector<uint64_t> detectors;
                std::set<uint64_t> frames;

                fp_t p = static_cast<fp_t>(inst.arg_data[0]);
                for (stim::DemTarget t : inst.target_data) {
                    if (t.is_relative_detector_id()) {
                        detectors.push_back(static_cast<uint64_t>(t.data + detector_offset));
                    } else if (t.is_observable_id()) {
                        frames.insert(static_cast<uint64_t>(t.data));
                    } else {
                        // This is just a separator, so call ef.
                        ef(p, detectors, frames);
                        detectors.clear();
                        frames.clear();
                    }
                }
                ef(p, detectors, frames);
            } else if (type == stim::DemInstructionType::DEM_SHIFT_DETECTORS) {
                detector_offset += static_cast<uint64_t>(inst.target_data[0].data);
            } else if (type == stim::DemInstructionType::DEM_DETECTOR) {
                for (stim::DemTarget t : inst.target_data) {
                    uint64_t d = static_cast<uint64_t>(t.data) + detector_offset;
                    df(d);
                }
            }
        }
    }
}

}   // graph
}   // qontra
