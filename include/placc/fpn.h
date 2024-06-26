/*
 *  author: Suhas Vittal
 *  date:   29 May 2024
 *
 *  Note that this is only for color codes.
 * */

#ifndef PLACC_FPN_h
#define PLACC_FPN_h

#include <qontra/graph.h>
#include <qontra/graph/algorithms/distance.h>
#include <qontra/graph/tanner_graph.h>

#include <qes.h>

namespace qg=qontra::graph;

namespace placc {

struct fpn_v_t : qg::base::vertex_t {
    enum class type { data, parity, flag };
    type qubit_type;
    bool is_widowed = false;
    // Parity qubit only data:
    std::vector<sptr<fpn_v_t>> support;
    std::map<std::pair<sptr<fpn_v_t>, sptr<fpn_v_t>>, sptr<fpn_v_t>>
        flag_usage_map;
};

struct fpn_e_t : qg::base::edge_t {};

typedef qg::DistanceMatrix<fpn_v_t, fp_t> DMAT;
typedef std::map<sptr<fpn_v_t>, size_t> mm_t;
typedef std::vector<std::vector<size_t>> ma_t;

// Class Stubs:
class ShorTree;

class FPN : public qg::Graph<fpn_v_t, fpn_e_t> {
public:
    FPN(qg::TannerGraph*, int color_to_remove);
    FPN(FPN&&) =default;

    void place_flags(void);
    void place_widowed_qubits(void);

    qes::Program<> phase_one_schedule(fp_t&, mm_t& x_ctr_map, mm_t& z_ctr_map);
    qes::Program<> phase_two_schedule(fp_t&, mm_t& x_ctr_map, mm_t& z_ctr_map, ma_t& x_flag_arr, ma_t& z_flag_arr);
private:
    void compute_cnot_order(void);

    DMAT compute_distance_matrix(const std::set<sptr<fpn_v_t>>& blocked_qubits);
    sptr<fpn_v_t> get_center_of(const std::vector<sptr<fpn_v_t>>& qubits,
                                const DMAT&,
                                const std::set<sptr<fpn_v_t>>& blocked_qubits);

    // Data structures for qubits and observables:
    std::vector<sptr<fpn_v_t>>  data_qubits,
                                inplace_parity_qubits,
                                removed_parity_qubits,
                                parity_qubits,
                                flag_qubits;
    std::vector<std::vector<sptr<fpn_v_t>>> obs_list;
    
    // Scheduling structures:
    // Entries: data/flag, time, parity
    vtils::TwoLevelMap<sptr<fpn_v_t>, size_t, sptr<fpn_v_t>> timestep_map;
    size_t max_timestep;

    uint64_t idctr;
};

void push_back_gate(qes::Program<>&, std::string, const std::vector<sptr<fpn_v_t>>&);
void push_back_gate(qes::Program<>&, std::string,
                        std::initializer_list<std::vector<sptr<fpn_v_t>>>);
void push_back_measurement(qes::Program<>&,
                            const std::vector<sptr<fpn_v_t>>>&, size_t& ctr, mm_t&);
void push_back_measurement(qes::Program<>&, 
                            std::initializer_list<std::vector<sptr<fpn_v_t>>>,
                            size_t& ctr, mm_t&);

}   // placc

#include "inl/fpn.inl"

#endif  // PLACC_FPN_h
