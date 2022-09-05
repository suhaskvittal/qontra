/*
 *  author: Suhas Vittal
 *  date:   27 August 2022
 * */

#ifndef GULLIVER_MEMORY_HIERARCHY_h
#define GULLIVER_MEMORY_HIERARCHY_h

#include "defs.h"
#include "decoding_graph.h"
#include "gulliver/defs.h"

#include <memory_system.h>

#include <chrono>
#include <map>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

/* We will implement a cache system ourselves.
 * Nothing too fancy, but gets the job done.
 * */

typedef uint64_t addr_t;

struct GulliverMemoryParams {
    uint n_sram_table_entries;
    uint n_detectors;

    std::string dram_config_file;
    std::string log_output_dir;
};

struct GulliverSramTableEntry {
    addr_t address;
    bool evictable;
};

class GulliverMemory {
public:
    GulliverMemory(const GulliverMemoryParams&);
    ~GulliverMemory();

    GulliverCycles access(uint, uint, bool mark_as_evictable);
    GulliverCycles prefetch(std::vector<uint>& detectors);

    void invalidate(void);

    uint64_t n_dram_accesses;
    uint64_t n_total_accesses;
private:
    uint64_t replace(addr_t);   // Strictly on chip.

    std::vector<GulliverSramTableEntry> sram_table;
    dramsim3::MemorySystem * main_memory;

    uint n_detectors;

    friend class Gulliver;
    friend void decoder_analysis_experiment(void);
};

uint
bound_detector(uint, uint n_detectors);
uint
unbound_detector(uint, uint n_detectors);

#endif
