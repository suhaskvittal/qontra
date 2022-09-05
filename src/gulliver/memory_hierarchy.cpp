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

GulliverMemory::GulliverMemory(const GulliverMemoryParams& params)
    :n_dram_accesses(0),
    n_total_accesses(0),
    sram_table(),
    main_memory(nullptr),
    n_detectors(params.n_detectors)
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

    GulliverSramTableEntry default_entry = {
        0x0,
        true
    };
    sram_table = std::vector<GulliverSramTableEntry>(params.n_sram_table_entries, 
                                                        default_entry);
}

GulliverMemory::~GulliverMemory() {
    delete main_memory;
}

GulliverCycles
GulliverMemory::access(uint di, uint dj, bool mark_as_evictable) {
    if (di > dj) {
        return access(dj, di, mark_as_evictable);
    }
    GulliverCycles n_cycles;

    uint bdi = bound_detector(di, n_detectors);
    uint bdj = bound_detector(dj, n_detectors);
    addr_t address = ADDRESS(bdi, bdj, n_detectors);
    // Examine sram table to check if entry exists.
    // If so, return it -- only one cycle lost!
    bool is_hit = false;
    for (GulliverSramTableEntry& e : sram_table) {
        if (e.address == address) {
            is_hit = true;
            if (mark_as_evictable) {
                e.evictable = true;
            }
            break;
        }
    }
    n_cycles.onchip++;
    n_total_accesses++;
    if (!is_hit) {
        n_dram_accesses++;
        // Then, we must go to DRAM.
        std::pair<uint, uint> di_dj = std::make_pair(di, dj);
        main_memory->AddTransaction(address, false); 
        while (memory_event_table.count(di_dj) == 0 ||
            !memory_event_table[di_dj])
        {
            main_memory->ClockTick(); 
            n_cycles.dram++;
        }
        memory_event_table[di_dj] = false;
        n_cycles.onchip += replace(address);
    }
    return n_cycles;
}

GulliverCycles
GulliverMemory::prefetch(std::vector<uint>& detectors) {
    invalidate();

    GulliverCycles n_cycles;
    
    std::map<std::pair<uint, uint>, addr_t> detector_pairs;
    for (uint ai = 0; ai < detectors.size(); ai++) {
        uint di = detectors[ai];
        uint bdi = bound_detector(di, n_detectors);
        for (uint bi = ai+1; bi < detectors.size(); bi++) {
            uint dj = detectors[bi];
            uint bdj = bound_detector(dj, n_detectors);
            // Compute address and tag
            addr_t address = ADDRESS(bdi, bdj, n_detectors);
            // First, check for table hit before prefetching
            // from memory.
            bool is_hit = false;
            for (GulliverSramTableEntry& e : sram_table) {
                if (e.address == address) {
                    e.evictable = false;
                    is_hit = true;
                    break;
                }
            }
            n_cycles.onchip++;
            if (!is_hit) {
                std::pair<uint, uint> di_dj = std::make_pair(di,dj);
                detector_pairs[di_dj] = address;
            }

            if (detector_pairs.size() == sram_table.size()) {
                goto prefetch_init_exit;
            }
        }
    }
prefetch_init_exit:
    // Now, we send memory requests to main memory when we can. We
    // wait until all requests are completed.
    for (auto di_dj_addr : detector_pairs) {
        auto di_dj = di_dj_addr.first;
        auto address = di_dj_addr.second;

        while (!main_memory->WillAcceptTransaction(address, false)) {
            main_memory->ClockTick();
            n_cycles.dram++;
        }
        main_memory->AddTransaction(address, false);
        main_memory->ClockTick();
        n_cycles.dram++;
        n_dram_accesses++;
        n_total_accesses++;
    }

    std::set<std::pair<uint, uint>> completed;
    do {
        for (auto di_dj_addr : detector_pairs) {
            auto di_dj = di_dj_addr.first;
            if (completed.count(di_dj)) {
                continue;
            }
            auto address = di_dj_addr.second;
            if (memory_event_table.count(di_dj) && memory_event_table[di_dj]) {
                replace(address);
                completed.insert(di_dj);
                memory_event_table[di_dj] = false;
            }
        }  
        main_memory->ClockTick();
        n_cycles.dram++;
    } while (completed.size() < detector_pairs.size());
    return n_cycles;
}

void
GulliverMemory::invalidate() {
    for (GulliverSramTableEntry& e : sram_table) {
        e.evictable = true;
    }
}

uint64_t
GulliverMemory::replace(addr_t address) {
    uint64_t n_cycles = 0;

    bool no_victim = true;
    uint victim;
    for (uint i = 0; i < sram_table.size(); i++) {
        GulliverSramTableEntry e = sram_table[i];
        if (e.evictable) {
            no_victim = false;
            victim = i;
        }
    }

    if (!no_victim) {
        // Spend a cycle evicting the entry.
        sram_table[victim] = (GulliverSramTableEntry) {
            address,
            false
        };
        n_cycles++;
    }
    
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
