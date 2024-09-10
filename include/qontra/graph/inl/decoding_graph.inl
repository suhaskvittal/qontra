/*
 *  author: Suhas Vittal
 *  date:   15 February 2024
 * */

namespace qontra {
namespace graph {

inline sptr<decoding::vertex_t>
DecodingGraph::get_boundary_vertex(int color) {
    uint64_t db = get_color_boundary_index(color);
    return get_vertex(db);
}

inline std::vector<sptr<decoding::hyperedge_t>>
DecodingGraph::get_all_edges() {
    return all_edges;
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
                size_t n_rep = static_cast<size_t>(inst.target_data[0].val());
                stim::DetectorErrorModel blk = dem.blocks[inst.target_data[1].val()];
                read_detector_error_model(blk, n_rep, detector_offset, ef, df);
            } else if (type == stim::DemInstructionType::DEM_ERROR) {
                std::vector<uint64_t> detectors;
                std::unordered_set<uint64_t> frames;

                fp_t p = static_cast<fp_t>(inst.arg_data[0]);
                for (stim::DemTarget t : inst.target_data) {
                    if (t.is_relative_detector_id()) {
                        detectors.push_back(t.val() + detector_offset);
                    } else if (t.is_observable_id()) {
                        frames.insert(t.val());
                    } else {
                        // This is just a separator, so call ef.
                        ef(p, detectors, frames);
                        detectors.clear();
                        frames.clear();
                    }
                }
                ef(p, detectors, frames);
            } else if (type == stim::DemInstructionType::DEM_SHIFT_DETECTORS) {
                detector_offset += inst.target_data[0].val();
            } else if (type == stim::DemInstructionType::DEM_DETECTOR) {
                for (stim::DemTarget t : inst.target_data) {
                    df(t.val() + detector_offset);
                }
            }
        }
    }
}

}   // graph
}   // qontra
