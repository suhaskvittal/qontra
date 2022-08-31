/*
 *  author: Suhas Vittal
 *  date:   27 August 2022
 * */

#include "gulliver/memory_hierarchy.h"

#define ADDRESS(di,dj,dn) (((di)*(dn)+(dj-di-1))*sizeof(fp_t))
#define GET_DI(x,dn) ((x) / (sizeof(fp_t)*dn))
#define GET_DJ(x,di,dn) ((x/sizeof(fp_t)) - ((di)*(dn)) + (di) + 1)

// Entry is true if it has been accessed.
static std::map<std::pair<uint, uint>, bool> memory_event_table;
// RNG for random replacement
static std::mt19937_64 CACHE_RNG(
        std::chrono::high_resolution_clock::now()
            .time_since_epoch().count());

GulliverCache::GulliverCache(const GulliverCacheParams& params)
    :n_accesses(0),
    n_misses(0),
    main_memory(nullptr),
    tag_store(),
    C(params.C),
    S(params.S),
    B(params.B),
    n_detectors(params.n_detectors),
    fake_cache(params.fake_cache)
{
    auto main_memory_callback = [this, memory_event_table] (uint64_t x) {
        uint di = GET_DI(x, this->n_detectors);
        uint dj = GET_DJ(x, di, this->n_detectors);
        di = unbound_detector(di, this->n_detectors);
        dj = unbound_detector(dj, this->n_detectors);
        std::pair<uint, uint> di_dj = std::make_pair(di, dj);
        memory_event_table[di_dj] = true;
    };
    // Finish initialization.
    main_memory = new dramsim3::MemorySystem(
            params.dram_config_file, params.log_output_dir,
            main_memory_callback, main_memory_callback);

    uint32_t n_sets = 1 << (C - B - S);
    uint32_t n_blks = 1 << S;
    GulliverCacheEntry default_entry = {
        0,
        0,
        0,
        0,
        0,
        false
    };
    tag_store = std::vector<std::vector<GulliverCacheEntry>>(n_sets, 
                    std::vector<GulliverCacheEntry>(n_blks, default_entry));
}

GulliverCache::~GulliverCache() {
    delete main_memory;
}

uint64_t
GulliverCache::access(uint di, uint dj) {
    if (di > dj) {
        return access(dj, di);
    }
    uint64_t n_cycles = 0;

    uint bdi = bound_detector(di, n_detectors);
    uint bdj = bound_detector(dj, n_detectors);
    addr_t address = ADDRESS(bdi, bdj, n_detectors);
    // Perform tag, index, and offset extraction (1 cycle).
    uint32_t offset = address & ((1 << B) - 1);
    uint32_t index = (address >> B) & ((1 << (C-B-S)) - 1);
    uint32_t tag = address >> (C-S);
    n_cycles++;
    // Access block in cache (1 cycle). 
    // If we miss, then we go to DRAM and evict data.
    // Note that we don't need to write to DRAM since
    // we don't modify any data in the cache.
    std::vector<GulliverCacheEntry> blocks = tag_store[index]; 
    bool is_hit = false;
    if (!fake_cache) {
        for (GulliverCacheEntry& blk : blocks) {
            if (blk.valid && blk.tag == tag) {
                is_hit = true;
                blk.last_use = 0;
                blk.n_accesses++;
            } else {
                blk.last_use++;
            }
            if (blk.is_new_entry) {
                blk.is_new_entry--;
            }
        }
        // Update stats.
        n_cycles++;
    }
    n_accesses++;
    if (!is_hit) {
        n_misses++;
        // Then, we must go to DRAM.
        std::pair<uint, uint> di_dj = std::make_pair(di, dj);
        main_memory->AddTransaction(address, false); 
        while (memory_event_table.count(di_dj) == 0 ||
            !memory_event_table[di_dj])
        {
            main_memory->ClockTick(); 
            n_cycles++;
        }
        memory_event_table[di_dj] = false;
        if (!fake_cache) {
            n_cycles += replace(address, tag, index, offset);
        }
    }
    return n_cycles;
}

uint64_t
GulliverCache::prefetch(std::vector<uint>& detectors) {
    if (fake_cache) {
        return 0;
    }

    // DEBUG: invalidate entire cache
//  for (auto& blks : tag_store) {
//      for (GulliverCacheEntry& blk : blks) {
//          blk.valid = false;
//      }
//  }
//  std::cout << "===========\n";

    uint64_t n_cycles = 0;

    std::map<std::pair<uint, uint>, addr_t> detector_pairs;
    for (uint ai = 0; ai < detectors.size(); ai++) {
        uint di = detectors[ai];
        uint bdi = bound_detector(di, n_detectors);
        for (uint bi = ai+1; bi < detectors.size(); bi++) {
            uint dj = detectors[bi];
            uint bdj = bound_detector(dj, n_detectors);
            // Compute address and tag
            addr_t address = ADDRESS(bdi, bdj, n_detectors);
            uint32_t offset = address & ((1 << B) - 1);
            uint32_t index = (address >> B) & ((1 << (C-B-S)) - 1);
            uint32_t tag = address >> (C-S);
//          std::cout << "[log] address(" << di << "," << dj << "): " 
//              << std::hex << address
//              << "\t [ " << tag << " | " 
//              << index << " | " << offset << " ]\n"
//              << std::dec;
            // First, check for cache hit before prefetching
            // from memory.
            std::vector<GulliverCacheEntry> blocks = tag_store[index]; 
            bool is_hit = false;
            for (GulliverCacheEntry& blk : blocks) {
                if (blk.valid && blk.tag == tag) {
                    is_hit = true;
                    blk.last_use = 0;
                    blk.n_accesses++;
                } else {
                    blk.last_use++;
                }
                if (blk.is_new_entry) {
                    blk.is_new_entry--;
                }
            }
            n_cycles++;
            if (!is_hit) {
                std::pair<uint, uint> di_dj = std::make_pair(di,dj);
                detector_pairs[di_dj] = address;
            }
        }
    }
    // Now, we send memory requests to main memory when we can. We
    // wait until all requests are completed.
    for (auto di_dj_addr : detector_pairs) {
        auto di_dj = di_dj_addr.first;
        auto address = di_dj_addr.second;
        uint32_t offset = address & ((1 << B) - 1);
        uint32_t index = (address >> B) & ((1 << (C-B-S)) - 1);
        uint32_t tag = address >> (C-S);

        while (!main_memory->WillAcceptTransaction(address, false)) {
            main_memory->ClockTick();
            n_cycles++;
        }
        main_memory->AddTransaction(address, false);
        main_memory->ClockTick();
        n_cycles++;
    }

    std::set<std::pair<uint, uint>> completed;
    std::set<uint32_t> new_tags;  // Track new tags so we don't replace 
                                  // too much.
    do {
        for (auto di_dj_addr : detector_pairs) {
            auto di_dj = di_dj_addr.first;
            if (completed.count(di_dj)) {
                continue;
            }
            auto address = di_dj_addr.second;
            if (memory_event_table.count(di_dj) && memory_event_table[di_dj]) {
                // Compute offset, index, tag for cache block.
                uint32_t offset = address & ((1 << B) - 1);
                uint32_t index = (address >> B) & ((1 << (C-B-S)) - 1);
                uint32_t tag = address >> (C-S);
                if (!fake_cache && new_tags.count(tag) == 0) {
                    replace(address, tag, index, offset);
                    new_tags.insert(tag);
                }
                completed.insert(di_dj);
                memory_event_table[di_dj] = false;
            }
        }  
        main_memory->ClockTick();
        n_cycles++;
    } while (completed.size() < detector_pairs.size());
    // And account for cycles for cache replacement:
    n_cycles += 2 * new_tags.size();
    return n_cycles;
}

uint64_t
GulliverCache::replace(addr_t address, uint64_t tag, uint64_t index,
        uint64_t offset) 
{
    uint64_t n_cycles = 0;
    // Assume victim selection takes one cycle.
    uint n_blocks = tag_store[index].size();

    uint victim;
    bool no_victim = true;
    uint32_t victim_last_use = 0;
    for (uint i = 0; i < n_blocks; i++) {
        GulliverCacheEntry blk = tag_store[index][i];
        if (!blk.valid || no_victim || blk.last_use > victim_last_use) {
            victim = i;
            no_victim = false;
            victim_last_use = blk.last_use;
        
            if (!blk.valid) {
                break;
            }
        }   
    }
    n_cycles++;
    // Assume replacement is one cycle.
    if (no_victim) {
        // Replace the first block.
        victim = 0;  
    }
//  if (tag_store[index][victim].valid && victim_last_use == 0) {
//      std::cout << "[error] couldn't find a man for the job?\n";
//  }
    tag_store[index][victim] = (GulliverCacheEntry) {
        address,
        tag,
        0,
        0,
        NEW_ENTRY_TTE,
        true
    };
    n_cycles++;
    return n_cycles;
}

uint
bound_detector(uint d, uint n_detectors) {
    return d == BOUNDARY_INDEX ? n_detectors-1 : d;
}

uint
unbound_detector(uint d, uint n_detectors) {
    return d == n_detectors-1 ? BOUNDARY_INDEX : d;
}
