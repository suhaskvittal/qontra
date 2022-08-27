/*
 *  author: Suhas Vittal
 *  date:   27 August 2022
 * */

#ifndef GULLIVER_MEMORY_HIERARCHY_h
#define GULLIVER_MEMORY_HIERARCHY_h

#include "defs.h"
#include "decoding_graph.h"

#include <memory_system.h>

#include <map>
#include <tuple>
#include <utility>
#include <vector>

/* We will implement a cache system ourselves.
 * Nothing too fancy, but gets the job done.
 * */

typedef uint64_t addr_t;

struct GulliverCacheParams {
    /* C, S, and B are standard cache parameters.
     * Facts:
     *  (1) 2**C = number of bytes in the cache.
     *  (2) S = 0 if the cache is direct mapped.
     *  (3) S = C - B if the cache is fully associative.
     *  (4) There are B bits in the block offset.
     *  (5) There are C - B - S in the index.
     *  (6) The remaining bits are given to the tag.
     * */
    uint C;
    uint S;
    uint B;
    uint n_detectors;

    uint tlbC;
    uint tlbB;

    std::string dram_config_file;
    std::string log_output_dir;
};

// Eviction grace period
#define NEW_ENTRY_TTE   2

struct GulliverCacheEntry {
    addr_t address;
    uint64_t tag;

    uint32_t n_accesses;
    uint8_t is_new_entry;
    bool valid;
};

/* We define a TLB-like cache for translating
 * detector index pairs into addresses. 
 * 
 * This will be fully associative and uses
 * the LRU policy. */
struct TLBEntry {
    uint di;
    uint dj;
    addr_t address;
    uint32_t last_use;
    bool valid;
};

class TLB {
public:
    TLB(uint C, uint B, uint n_detectors);
    /* First part is address, second part is
     * a boolean which is true if there was a miss,
     * and third part is the number of cycles taken.*/
    typedef std::tuple<addr_t, bool, uint64_t> AddressResult;
    AddressResult address(uint, uint);
private:
    uint64_t replace(addr_t, uint, uint);
    /* No need for set-level access. */
    std::vector<TLBEntry> tag_store;

    uint C;
    uint B;
    uint n_detectors;
};


class GulliverCache {
public:
    GulliverCache(const GulliverCacheParams&);
    ~GulliverCache();

    uint64_t access(uint, uint);

    // Statistics
    uint32_t n_accesses;
    uint32_t n_misses;
    uint32_t n_tlb_accesses;
    uint32_t n_tlb_misses;
private:
    uint64_t replace(addr_t, uint64_t tag, uint64_t index, uint64_t offset);
    dramsim3::MemorySystem * main_memory;
    TLB * tlb;

    /* Tag store is accessed by:
     * (1) Set
     * (2) Way
     * */
    std::vector<std::vector<GulliverCacheEntry>> tag_store;

    uint C;
    uint S;
    uint B;
};

uint
bound_detector(uint, uint n_detectors);
uint
unbound_detector(uint, uint n_detectors);

#endif
