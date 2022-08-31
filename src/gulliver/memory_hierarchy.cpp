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
    :main_memory(nullptr),
    sram_table(),
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
        false
    };
    sram_table = std::vector<GulliverSramTableEntry>(params.n_sram_table_entries, 
                                                        default_entry);
}

GulliverMemory::~GulliverMemory() {
    delete main_memory;
}

uint64_t
GulliverMemory::access(uint di, uint dj, bool mark_as_evictable) {
    if (di > dj) {
        return access(dj, di, mark_as_evictable);
    }
    uint64_t n_cycles = 0;

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
    n_cycles++;
    if (!is_hit) {
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
        n_cycles += replace(address);
    }
    return n_cycles;
}

uint64_t
GulliverMemory::prefetch(std::vector<uint>& detectors) {
    invalidate();

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
            // First, check for table hit before prefetching
            // from memory.
            bool is_hit = false;
            for (GulliverSramTableEntry& e : sram_table) {
                if (e.address == address) {
                    is_hit = true;
                    break;
                }
            }
            n_cycles++;
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
            n_cycles++;
        }
        main_memory->AddTransaction(address, false);
        main_memory->ClockTick();
        n_cycles++;
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
        n_cycles++;
    } while (completed.size() < detector_pairs.size());
    return n_cycles;
}

uint64_t
GulliverMemory::invalidate() {
    uint64_t n_cycles = 0;
    for (GulliverSramTableEntry& e : sram_table) {
        e.evictable = true;
    }
    n_cycles++;
    return n_cycles;
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
