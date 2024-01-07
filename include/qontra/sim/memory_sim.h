/*
 *  author: Suhas Vittal
 *  date:   10 December 2023
 * */

#ifndef QONTRA_SIM_MEMORY_h
#define QONTRA_SIM_MEMORY_h

#include "qontra/graph/lattice_graph.h"
#include "qontra/sim/frame_sim.h"
#include "qontra/stats.h"
#include "qontra/tables.h"

#include <map>
#include <string>
#include <vector>

namespace qontra {

class MemorySimulator {
public:
    MemorySimulator(graph::LatticeGraph&);

    void    initialize(void);
    void    reset(void);

    void    run(uint64_t shots);

    enum class lrc_policy_t { none, always, optimal, eraser };
    enum class lrc_circuit_t { swap, dqlr };
    struct {
        uint distance;
        uint rounds;
        ErrorTable  errors;
        TimeTable   timing;

        std::string stim_output_file;
        std::string syndrome_output_folder;
        std::string data_output_file;

        bool is_memory_x=false;
        bool enable_leakage=false;
        
        lrc_policy_t    lrc_policy = lrc_policy_t::none;
        lrc_circuit_t   lrc_circuit = lrc_circuit_t::swap;
        size_t          lrc_stride_size = 1;
    } config;

    // Statistics:
    statistic_t<fp_t>   s_lrcs_per_round;
    statistic_t<fp_t>   s_leakage_population_ratio;
private:
    typedef fp_t time_t;

    void    run_batch(uint64_t shots_in_batch);

    // do_gate and do_measurement accept an optional trial argument. By default, trial = -1,
    // which means both functions perform the operation across all trials and inject errors
    // stochastically.
    //
    // When trial >= 0, then we will only perform an operation (and the corresponding errors)
    // for that specified trial.
    time_t  do_gate(std::string op, std::vector<uint> operands, int64_t trial=-1);
    time_t  do_measurement(std::vector<uint> operands, int64_t trial=-1);
    void    inject_idling_error_positive(std::vector<uint> on_qubits, int64_t trial=-1);
    void    inject_idling_error_negative(std::vector<uint> not_on_qubits, int64_t trial=-1);

    // The below functions are deterministic and should not vary from trial to trial.
    // But, they may leverage trial-specific information (i.e. timing error may be different
    // between trials).
    void    create_event_or_obs(std::vector<uint> operands, bool create_event=true);
    void    inject_timing_error(std::vector<uint> qubits);
    
    graph::LatticeGraph lattice_graph;
    uint64_t n_qubits;
    uint64_t n_detection_events;
    uint64_t n_observables;

    std::vector<uint>   all_qubits;
    std::vector<uint>   data_qubits;
    std::vector<uint>   parity_qubits;
    std::vector<uint>   xp_qubits;
    std::vector<uint>   zp_qubits;

    std::vector<uint>   syndrome_parity_qubits;

    std::vector<std::vector<uint>>  obs_list;

    uptr<FrameSimulator> sim;
    stim::simd_bit_table<SIMD_WIDTH> syndromes;
    // Tracking data structures.
    time_t  elapsed_time;
    std::map<uint64_t, time_t> shot_time_delta_map;

    uint64_t meas_ctr;
    uint64_t event_ctr;
    uint64_t obs_ctr;
    std::map<uint, uint64_t> meas_ctr_map;
    std::map<uint, uint64_t> event_ctr_map;

    stim::Circuit   sample_circuit;
    bool            is_recording_stim_instructions;

    // 
    // LRC functions:
    //
    void    lrc_reset(void);
    void    lrc_execute_lrcs_from_await_queue(void);

    stim::simd_bits<SIMD_WIDTH> lrc_optimal_identify_lrcs(void);
                                    // Returns one for shots where an LRC should be used.
    void                        lrc_optimal_perform_lrcs(void); 

    void    lrc_measure_qubits(const std::map<uint, uint>& swap_set, int64_t trial=-1);

    std::map<uint, uint>    lrc_solve_maximum_matching(
                                    const std::vector<uint>& avail_data,
                                    const std::vector<uint>& avail_parity);
    //
    // LRC variables:
    //
    std::deque<uint> lrc_await_queue;   // Only used by "always" LRCs.
    std::map<uint64_t, std::map<uint, uint>> lrc_optimal_lrc_map_table;
    //
    // ERASER functions:
    //
    void    eraser_initialize(void);
    void    eraser_reset(void);

    stim::simd_bits<SIMD_WIDTH> eraser_execute_lrcs(void);   
                                    // Returns 1 wherever an LRC was used.
    void    eraser_update_syndrome_buffer(const std::map<uint, uint64_t>& prev_meas_ctr_map);
    //
    // ERASER variables:
    //
    std::map<uint, uint> eraser_qubit_to_syndrome_buffer_index_map;
    std::map<uint, std::array<uint, 2>> eraser_swap_lookup_table;
    // The recently_scheduled_qubits_table should be #shots by #qubits for better cache usage.
    stim::simd_bit_table<SIMD_WIDTH>    eraser_recently_scheduled_qubits_table;
    stim::simd_bit_table<SIMD_WIDTH>    eraser_syndrome_buffer;
};

}   // qontra

#endif  // QONTRA_SIM_MEMORY_h
