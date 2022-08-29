/*
 *  author: Suhas Vittal
 *  date:   27 August 2022
 * */

#include "gulliver/memory_hierarchy.h"

// Entry is true if it has been accessed.
static std::map<std::pair<uint, uint>, bool> memory_event_table;
// RNG for random replacement
static std::mt19937_64 CACHE_RNG(
        std::chrono::high_resolution_clock::now()
            .time_since_epoch().count());

GulliverCache::GulliverCache(const GulliverCacheParams& params)
    :n_accesses(0),
    n_misses(0),
    n_tlb_accesses(0),
    n_tlb_misses(0),
    main_memory(nullptr),
    tlb(nullptr),
    tag_store(),
    C(params.C),
    S(params.S),
    B(params.B),
    fake_cache(params.fake_cache)
{
    uint n_detectors = params.n_detectors;
    auto main_memory_callback = [n_detectors, memory_event_table] (uint64_t x) {
        uint di = x / n_detectors;
        uint dj = x - di * n_detectors;
        di = unbound_detector(di, n_detectors);
        dj = unbound_detector(dj, n_detectors);
        std::pair<uint, uint> di_dj = std::make_pair(di, dj);
        memory_event_table[di_dj] = true;
    };
    // Finish initialization.
    main_memory = new dramsim3::MemorySystem(
            params.dram_config_file, params.log_output_dir,
            main_memory_callback, main_memory_callback);
    tlb = new TLB(params.tlbC, params.tlbB, n_detectors);

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
    delete tlb;
}

uint64_t
GulliverCache::access(uint di, uint dj) {
    uint64_t n_cycles = 0;

    TLB::AddressResult tlb_result = tlb->address(di, dj);
    addr_t address = std::get<0>(tlb_result);
    // Update stats.
    n_tlb_accesses++;
    n_tlb_misses += std::get<1>(tlb_result) ? 0 : 1;
    n_cycles += std::get<2>(tlb_result);
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
                if (blk.is_new_entry) {
                    blk.is_new_entry--;
                }
                blk.n_accesses++;
            } else {
                blk.last_use++;
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
GulliverCache::replace(addr_t address, uint64_t tag, uint64_t index,
        uint64_t offset) 
{
    uint64_t n_cycles = 0;
    // Assume victim selection takes one cycle.
    uint n_blocks = tag_store[index].size();
    bool no_victim = true;
    uint victim;
#ifdef GC_POLICY_LRU
    uint32_t victim_last_use = 0;
    for (uint i = 0; i < n_blocks; i++) {
        GulliverCacheEntry blk = tag_store[index][i];
        if (!blk.valid || no_victim || blk.last_use > victim_last_use) {
            victim = i;
            no_victim = false;
            victim_last_use = blk.last_use;
        }   
    }
#elifdef GC_POLICY_RR
    auto bi = CACHE_RNG() % n_blocks;
    auto bj = CACHE_RNG() % n_blocks;
    GulliverCacheEntry blk_i = tag_store[index][bi];
    GulliverCacheEntry blk_j = tag_store[index][bj];
    if (blk_i.last_use > blk_j.last_use) {
        victim = bi;
    } else {
        victim = bj;
    }
    no_victim = false; 
#else
    uint32_t victim_accesses = 0;
    for (uint i = 0; i < n_blocks; i++) {
        GulliverCacheEntry blk = tag_store[index][i];
        if (blk.is_new_entry) {
            continue;
        }
        if (!blk.valid || no_victim || blk.n_accesses < victim_accesses) {
            victim = i;
            no_victim = false;
            victim_accesses = blk.n_accesses;
        }
    }
#endif
    n_cycles++;
    // Assume replacement is one cycle.
    if (no_victim) {
        // Replace the first block.
        victim = 0;  
    }
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

TLB::TLB(uint C, uint B, uint n_detectors)
    :tag_store((1<<(C-B)), (TLBEntry){0,0,0,0,false}),
    C(C),
    B(B),
    n_detectors(n_detectors)
{}

TLB::AddressResult
TLB::address(uint di, uint dj) {
    if (di > dj) {
        // Enforce di <= dj.
        return address(dj, di);
    }
    uint64_t n_cycles = 0;
    // TLB access is 1 cycle.
    bool is_hit = false;
    addr_t address = 0x0;
    for (TLBEntry& blk : tag_store) {
        if (blk.valid && blk.di == di && blk.dj == dj) {
            is_hit = true; 
            address = blk.address;
            blk.last_use = 0;
        } else {
            blk.last_use++;
        }
    }
    n_cycles++;
    if (!is_hit) {
        // Then, we spend three cycles (conservatively)
        // on "address translation".
        
        // Bound detectors to help with address translation.
        uint bdi = bound_detector(di, n_detectors);
        uint bdj = bound_detector(dj, n_detectors);
        // Multiplication is obviously not one cycle,
        // but we can optimize it in hardware because 
        // n_detectors would be fixed. We could take the number
        // of bits in the power of 2 larger than n_detectors.
        // and do a bit shift.
        address = bdi * n_detectors + bdj;
        n_cycles += 3;
        // And a few more on replacement.
        n_cycles += replace(address, di, dj);
    }
    return std::make_tuple(address, is_hit, n_cycles);
}

uint64_t
TLB::replace(addr_t address, uint di, uint dj) {
    // Identifying the victim is one cycle.
    uint64_t n_cycles = 0; 
    uint n_blocks = tag_store.size();
    uint victim = 0;
    uint32_t victim_last_use = tag_store[0].last_use;
    for (uint i = 1; i < tag_store.size(); i++) {
        TLBEntry blk = tag_store[i];
        if (!blk.valid || blk.last_use > victim_last_use) {
            victim = i;
            victim_last_use = blk.last_use;
        }
    }
    n_cycles++;
    // Eviction + installation is one cycle.
    tag_store[victim] = (TLBEntry) {
        di,
        dj,
        address,
        0,
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
