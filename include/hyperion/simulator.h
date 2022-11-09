/*
 *  author: Suhas Vittal
 *  date:   7 September 2022
 * */

#ifndef HYPERION_SIMULATOR_h
#define HYPERION_SIMULATOR_h

#include "defs.h"
#include "decoding_graph.h"
#include "mwpm_decoder.h"

#include <memory_system.h>

#include <algorithm>
#include <array>
#include <deque>
#include <map>
#include <vector>
#include <utility>

//#define HSIM_DEBUG
#define FILTER_CUTOFF 10
#define WINDOW_SIZE 100

namespace qrc {
namespace hyperion {

struct HyperionSimulatorParams {
    uint n_detectors;
    uint n_detectors_per_round;

    uint n_registers;
    uint bfu_fetch_width;
    uint bfu_compute_stages;
    uint bfu_priority_queue_size;

    uint8_t bankgroup;
    uint8_t bank;
    uint32_t row_offset;
};

struct MemoryEventEntry {
    uint di;
    uint dj;
    uint logical_qubit;
};

class HyperionDeque;

class HyperionSimulator {
public:
    HyperionSimulator(dramsim3::MemorySystem*,
            std::map<addr_t, bool> * memory_event_table,
            const PathTable& path_table,
            const HyperionSimulatorParams&);

    bool load_detectors(const std::vector<uint>&);  // Returns false if the
                                                    // syndrome is not 
                                                    // servicable.
    void load_path_table(const PathTable&);
    void load_qubit_number(uint);
    void load_base_address(uint8_t bankgroup, uint8_t bank, uint32_t row_offset);

    void tick(void);
    void sig_end_round(uint=1);  // Number of rounds to jump. 
    
    bool is_idle(void);
    void force_idle(void);
    std::map<uint, uint> get_matching(void);

    void reset_stats(void);
    
    dramsim3::MemorySystem * dram;
    // Statistics
    uint64_t rowhammer_flips(void);
    uint64_t row_activations(void);

    uint64_t prefetch_cycles;
    uint64_t bfu_cycles;
protected:
    void tick_prefetch(void);
    void tick_bfu(void);

    void tick_bfu_compute(uint stage);
    void tick_bfu_sort(uint stage);
    void tick_bfu_fetch(void);

    bool access(addr_t, bool set_evictable_on_hit);
    void update_state(void);

    void clear(void);

    enum class State { prefetch, bfu, idle };

    struct Register {
        addr_t address;
        uint64_t last_use;
        bool valid;
    };

    struct DequeEntry {
        std::map<uint, uint> running_matching;
        uint32_t matching_weight;
        uint next_unmatched_index;
        
        bool operator<(const DequeEntry& other) const {
            return matching_weight/running_matching.size() 
                < other.matching_weight/other.running_matching.size();
        }
    };

    struct BFUPipelineLatch {
        std::deque<std::pair<uint, uint32_t>> proposed_matches;
        DequeEntry base_entry;
        bool stalled;
        bool valid;
    };

/* Microarchitectural components.*/
    // Global memory
    std::vector<Register> register_file;    
    std::vector<addr_t> dram_await_array;
    std::vector<uint> detector_vector_register; // Holds the detectors from the
                                                // current syndrome.
    uint32_t mean_weight_register;
    uint access_counter;
// Prefetch
    uint min_detector_register;
    uint major_detector_register;
    std::map<uint, uint> minor_detector_table;
    // There are fetch_width hardware deques of size 2*fetch_width
    std::vector<HyperionDeque> hardware_deques;
    // Size of latches is fetch_width by (1 + 4 + 1)
    // There is one FETCH stage,
    //          four SORT stages,
    //          and 1 compute stage.
    std::vector<std::vector<BFUPipelineLatch>> bfu_pipeline_latches;   
    // Replacement policy
    DequeEntry best_matching_register;
    std::deque<addr_t> replacement_queue;
    // Global state machine
    State state; 
    bool bfu_idle;

    /* Data */
    std::map<addr_t, bool> * memory_event_table;
    PathTable path_table;
    /* Configuration parameters. */
private:
    uint curr_max_detector;

    uint n_detectors;
    uint n_detectors_per_round;

    uint bfu_fetch_width;
    uint bfu_compute_stages; 

    uint curr_qubit;
    addr_t base_address;

    bool has_boundary;

    friend class HyperionDeque;
};

// Fixed-size priority queue.
class HyperionDeque {
public:
    HyperionDeque(uint);

    void push(HyperionSimulator::DequeEntry);
    HyperionSimulator::DequeEntry top(void);
    void pop(void);
    uint size(void);
    bool empty(void);
    void clear(void);
private:
    void downheap(uint);
    void upheap(uint);

    uint parent(uint);
    uint left(uint);
    uint right(uint);

    std::vector<HyperionSimulator::DequeEntry> backing_array;

    uint max_size;

    friend class HyperionSimulator;
};

addr_t get_base_address(uint8_t bankgroup, uint8_t bank, uint32_t row_offset,
        dramsim3::Config*);
addr_t to_address(uint, uint, addr_t base, uint n_detectors);
std::pair<uint, uint> from_address(addr_t, addr_t base, uint n_detectors);

uint bound_detector(uint, uint n_detectors);
uint unbound_detector(uint, uint n_detectors);

} // hyperion
} // qrc

#endif
