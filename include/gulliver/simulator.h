/*
 *  author: Suhas Vittal
 *  date:   7 September 2022
 * */

#ifndef GULLIVER_SIMULATOR_h
#define GULLIVER_SIMULATOR_h

#include "defs.h"
#include "decoding_graph.h"

#include <memory_system.h>

#include <deque>
#include <map>
#include <stack>
#include <vector>
#include <utility>

struct GulliverSimulatorParams {
    uint n_detectors;

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
            const std::map<std::pair<uint, uint>, fp_t>& weight_table,
            const GulliverSimulatorParams&);

    void load(const std::vector<uint>&);
    void tick(void);
    
    bool is_idle(void);
    std::map<uint, uint> get_matching(void);
    
    // Statistics
    uint64_t rowhammer_flips(void);
    uint64_t row_activations(void);

    uint64_t prefetch_cycles;
    uint64_t bfu_cycles;
protected:
    void tick_prefetch(void);
    void tick_predecode(void);
    void tick_bfu(void);

    void tick_bfu_compute(uint stage);
    void tick_bfu_fetch();

    void clear(void);

    enum State { prefetch, predecode, bfu, idle };

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

    struct PipelineLatch {
        std::stack<std::pair<uint, fp_t>> proposed_matches;
        StackEntry base_entry;
        bool stalled;
        bool valid;
    };

    /* Microarchitectural components.*/
    dramsim3::MemorySystem * dram;
    std::vector<Register> register_file;    
    std::pair<uint, uint> next_dram_request_register;
    std::vector<std::pair<uint, uint>> dram_await_array;

    std::vector<uint> detector_vector_register; // Holds the detectors from the
                                                // current syndrome.
    std::stack<StackEntry> hardware_stack;
    std::vector<std::vector<PipelineLatch>> pipeline_latches;   // Size is 
                                                                // fetch_width by
                                                                // (hw_threshold-1)
    StackEntry best_matching_register;
    std::deque<addr_t> replacement_queue;
        
    State state; 
    bool bfu_idle;

    /* Data */
    std::map<addr_t, bool> * memory_event_table;

    std::map<std::pair<uint, uint>, fp_t> weight_table;
    /* Configuration parameters. */
private:
    uint n_detectors;

    uint bfu_fetch_width;
    uint bfu_hw_threshold; 

    uint8_t bankgroup;
    uint8_t bank;
    uint32_t row_offset;

    addr_t base_address;
};

addr_t get_base_address(uint8_t bankgroup, uint8_t bank, uint32_t row_offset,
        dramsim3::Config*);
addr_t to_address(uint, uint, addr_t base, uint n_detectors);
std::pair<uint, uint> from_address(addr_t, addr_t base, uint n_detectors);

uint bound_detector(uint, uint n_detectors);
uint unbound_detector(uint, uint n_detectors);

#endif
