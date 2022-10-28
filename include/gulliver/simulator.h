/*
 *  author: Suhas Vittal
 *  date:   27 October 2022
 * */

#ifndef GULLIVER_SIMULATOR_h
#define GULLIVER_SIMULATOR_h

#include "defs.h"
#include "decoding_graph.h"
#include "mwpm_decoder.h"

#include <memory_system.h>

#include <algorithm>
#include <deque>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace qrc {
namespace gulliver {

#define GSIM_DEBUG

class GulliverSimulator {
public:
    GulliverSimulator(dramsim3::MemorySystem*,
            std::map<addr_t, bool> * memory_event_table,
            const DecodingGraph&,
            MWPMDecoder*,
            uint search_width);

    void load_syndrome(const std::vector<uint8_t>&);
    void tick(void);
    bool is_idle(void);
    void force_idle(void);
    std::map<uint, uint> get_matching(void);

    void reset_stats(void);
    
    dramsim3::MemorySystem * dram;

    uint64_t traverse_cycles;
    uint64_t combine_cycles;
    uint64_t max_active_syndromes_in_box;
protected:
    void tick_traverse(void);
    void tick_combine(void);

    void tick_traverse_mem(void);
    void tick_traverse_read(void);

    void update_state(void);
    void clear(void);

    enum class State {traverse, combine, idle}; 

    struct TraversePipelineLatch {
        std::vector<uint> detector_array;
        std::map<uint, uint> matching;
        addr_t dram_address;
        bool valid;
        bool is_stalled;
    };

    std::vector<uint8_t> syndrome_register;
    uint hamming_weight_register;
    std::vector<addr_t> dram_await_array;
    State state;

    // TRAVERSE
    uint current_qubit;     // Lower left corner of box.
    std::vector<TraversePipelineLatch> traverse_pipeline_latches;
    std::map<uint, uint> cumulative_matching_register;
    std::map<uint, std::set<uint>> matchings;
    bool traverse_idle;

    // Non-microarchitectural structures (to help simulation)
    std::map<addr_t, std::vector<uint8_t>> addr_to_syndrome;  // Do MWPM after DRAM
                                                              // access finishes.
    MWPMDecoder * decoder;  // Note that the decoder is not used
                            // for the true design. We are using it
                            // to simulate the results from memory
                            // accesses to DRAM.
    DecodingGraph decoding_graph;
private:
    std::map<addr_t, bool> * memory_event_table;

    uint search_width;
};

} //gulliver
} // qrc

#endif
