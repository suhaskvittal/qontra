/*
 *  author: Suhas Vittal
 *  date:   7 September 2022
 * */

#ifndef GULLIVER_SIMULATOR_h
#define GULLIVER_SIMULATOR_h

#include "defs.h"
#include "decoding_graph.h"

#include <memory_system.h>

#include <boost/graph/graph_traits.hpp>

#include <deque>
#include <map>
#include <stack>
#include <vector>
#include <utility>

#define N_ROWS_ADJACENCY_LISTS 8
#define FILTER_CUTOFF 8

struct GulliverSimulatorParams {
    uint n_detectors;
    uint n_detectors_per_round;

    uint n_registers;
    uint bfu_fetch_width;
    uint bfu_hw_threshold;

    uint8_t bankgroup;
    uint8_t bank;
    uint32_t row_offset;
};

struct MemoryEventEntry {
    uint di;
    uint dj;
    uint logical_qubit;
};

class GulliverSimulator {
public:
    GulliverSimulator(dramsim3::MemorySystem*, 
            std::map<addr_t, bool> * memory_event_table,
            const PathTable& path_table,
            const GulliverSimulatorParams&);

    bool load_detectors(const std::vector<uint>&);  // Returns false if the
                                                    // syndrome is not 
                                                    // servicable.
    void load_path_table(const PathTable&);
    void load_base_address(uint8_t bankgroup, uint8_t bank, uint32_t row_offset);

    void tick(void);
    
    bool is_idle(void);
    std::map<uint, uint> get_matching(void);

    void reset_stats(void);
    
    // Statistics
    uint64_t rowhammer_flips(void);
    uint64_t row_activations(void);

    uint64_t prefetch_cycles;
    uint64_t bfu_cycles;
protected:
    void setup(void);

    void tick_prefetch(void);
    void tick_bfu(void);

    void tick_bfu_compute(uint stage);
    void tick_bfu_fetch();

    bool access(addr_t, bool set_evictable_on_hit);

    void clear(void);

    enum State { prefetch, bfu, idle };

    struct Register {
        addr_t address;
        bool evictable;
        bool valid;
    };

    struct StackEntry {
        std::map<uint, uint> running_matching;
        fp_t matching_weight;
        uint next_unmatched_index;
    };

    struct PDScoreboardEntry {
        uint best_mate_index;
        uint n_suitors;
        fp_t mate_weight;
    };

    struct BFUPipelineLatch {
        std::stack<std::pair<uint, fp_t>> proposed_matches;
        StackEntry base_entry;
        bool stalled;
        bool valid;
    };

    /* Microarchitectural components.*/
    // Global memory
    dramsim3::MemorySystem * dram;
    std::vector<Register> register_file;    
    std::pair<uint, uint> next_dram_request_register;
    std::vector<addr_t> dram_await_array;
    std::vector<uint> detector_vector_register; // Holds the detectors from the
                                                // current syndrome.
    // Predecode
    std::vector<PDScoreboardEntry> predecode_scoreboard;
    uint index_register;
    uint neighbor_register;
    // BFU
    std::stack<StackEntry> hardware_stack;
    // Size of latches is fetch_width by (hw_threshold-1)
    std::vector<std::vector<BFUPipelineLatch>> bfu_pipeline_latches;   
    // Replacement policy
    StackEntry best_matching_register;
    std::deque<addr_t> replacement_queue;
    // Global state machine
    State state; 
    bool bfu_idle;

    /* Data */
    std::map<addr_t, bool> * memory_event_table;
    PathTable path_table;
    /* Configuration parameters. */
private:
    uint n_detectors;
    uint n_detectors_per_round;

    uint critical_index;

    uint bfu_fetch_width;
    uint bfu_hw_threshold; 

    addr_t base_address;
};

addr_t get_base_address(uint8_t bankgroup, uint8_t bank, uint32_t row_offset,
        dramsim3::Config*);
addr_t to_address(uint, uint, addr_t base, uint n_detectors);
std::pair<uint, uint> from_address(addr_t, addr_t base, uint n_detectors);

uint bound_detector(uint, uint n_detectors);
uint unbound_detector(uint, uint n_detectors);

#endif
