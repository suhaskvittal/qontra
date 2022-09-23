/*
 *  author: Suhas Vittal
 *  date:   7 September 2022
 * */

#include "gulliver/simulator.h"
#include <limits>

//#define GSIM_DEBUG

namespace qrc {
namespace gulliver {

GulliverSimulator::GulliverSimulator(dramsim3::MemorySystem * dram, 
        QubitCache * cache,
        std::map<addr_t, bool> * memory_event_table,
        const PathTable& path_table,
        const GulliverSimulatorParams& params)
    :
    /* Statistics */
    prefetch_cycles(0),
    bfu_cycles(0),
    /* Microarchitecture */
    dram(dram),
    cache(cache),
    register_file(params.n_registers, (Register){0x0, 0, false}),
    dram_await_array(),
    detector_vector_register(),
    ccomp_detector_register(0),
    major_detector_register(0),
    minor_detector_table(),
    hardware_stack(),
    bfu_pipeline_latches(params.bfu_fetch_width, 
            std::vector<BFUPipelineLatch>(params.bfu_hw_threshold-1)),
    best_matching_register(),
    replacement_queue(),
    state(GulliverSimulator::State::prefetch),
    bfu_idle(false),
    /* Data */
    memory_event_table(memory_event_table),
    path_table(path_table),
    /* Config parameters */
    curr_max_detector(0),
    n_detectors(params.n_detectors),
    n_detectors_per_round(params.n_detectors_per_round),
    bfu_fetch_width(params.bfu_fetch_width),
    bfu_hw_threshold(params.bfu_hw_threshold),
    curr_qubit(0),
    base_address(0)
{
    clear();
    load_base_address(params.bankgroup, params.bank, params.row_offset);
}

bool
GulliverSimulator::load_detectors(const std::vector<uint>& detector_array) {
    clear();
    detector_vector_register = detector_array;
    // Compute the critical index.
    curr_max_detector = n_detectors_per_round;

    bool computable = detector_array.size() <= bfu_hw_threshold;
    if (!computable) {
        state = State::idle;  // Don't even bother trying.
    }
    return computable;
}

void
GulliverSimulator::load_path_table(const PathTable& pt) {
    path_table = pt;  
}

void
GulliverSimulator::load_qubit_number(uint qubit) {
    curr_qubit = qubit;
}

void
GulliverSimulator::load_base_address(uint8_t bg, uint8_t ba, uint32_t ro_offset) {
    base_address = get_base_address(bg, ba, ro_offset, dram->config_);
}

void
GulliverSimulator::tick() {
    if (cache != nullptr) {
        cache->tick();
    }
    // Update last use for every register.
    for (Register& r : register_file) {
        r.last_use++;
    }
    // Go perform any replacements if necessary.
    if (replacement_queue.size()) {
        addr_t address = replacement_queue.front();
        // Check if there are any evictable registers.
        Register * victim = &register_file[0];
        for (uint i = 1; i < register_file.size(); i++) {
            Register& r = register_file[i];
            if (victim->valid &&
                    (!r.valid || r.last_use > victim->last_use)) 
            {
                victim = &register_file[i];
            }
        }
        victim->address = address;
        victim->valid = true;
        victim->last_use = 0;

        replacement_queue.pop_front();
#ifdef GSIM_DEBUG
        std::cout << "\t[Gulliver] Installed " << std::hex 
            << address << std::dec << "\n";
#endif
    }
    // Check if any transactions to DRAM have completed.
    for (auto it = dram_await_array.begin(); it != dram_await_array.end(); ) {
        addr_t address = *it;
        if (memory_event_table->count(address) 
                && memory_event_table->at(address)) 
        {
            it = dram_await_array.erase(it);
            memory_event_table->at(address) = false;
            // Add to replacement queue.
            auto di_dj = from_address(address, base_address, n_detectors);
#ifdef GSIM_DEBUG
            std::cout << "\t[DRAM] Retrieved " << std::hex
                << address << std::dec << "(" << di_dj.first
                << "," << di_dj.second << ")\n";
#endif
            uint di = di_dj.first;
            if (state != State::ccomp) {
                replacement_queue.push_back(address);
            }
        } else {
            it++;
        }
    }

    switch (state) {
    case State::ccomp:
        ccomp_cycles++;
        tick_ccomp();
        break;
    case State::prefetch:
        prefetch_cycles++;
        tick_prefetch();
        break;
    case State::bfu:
        bfu_cycles++;
        tick_bfu();
        break;
    default:
        break;
    }
}

void
GulliverSimulator::sig_end_round() {
    // We subtract by 1 because we don't count the boundary.
    if (curr_max_detector < n_detectors-1) {
        curr_max_detector += n_detectors_per_round;
    }
#ifdef GSIM_DEBUG
    std::cout << "SIGNAL -- ROUND END\n";
#endif
}

bool
GulliverSimulator::is_idle() {
    return state == GulliverSimulator::State::idle;
}

std::map<uint, uint>
GulliverSimulator::get_matching() {
    return best_matching_register.running_matching;
}

void
GulliverSimulator::reset_stats() {
    prefetch_cycles = 0;
    bfu_cycles = 0;
}

uint64_t
GulliverSimulator::rowhammer_flips() {
    return dram->dram_system_->rowhammer_flips();
}

uint64_t
GulliverSimulator::row_activations() {
    return dram->dram_system_->row_activations();
}

void 
GulliverSimulator::tick_ccomp() {
    if (cache->completion_queue.empty()) {
        update_state();
        return;
    }
    // Retrieve data from memory.
    uint di = cache->completion_queue.front();
    uint dj = unbound_detector(ccomp_detector_register, n_detectors);
    if (di == dj) {
        // Update ccomp detector register
        ccomp_detector_register++;
        if (ccomp_detector_register >= n_detectors) {
            ccomp_detector_register = 0;
            cache->complete(curr_qubit, di);
            cache->completion_queue.pop_front();
        }
    }
    addr_t address = to_address(di, dj, base_address, n_detectors);
    if (dram->WillAcceptTransaction(address, false)) {
        dram->AddTransaction(address, false);
        dram_await_array.push_back(address);
        // Update ccomp detector register
        ccomp_detector_register++;
        if (ccomp_detector_register >= n_detectors) {
            ccomp_detector_register = 0;
            cache->complete(curr_qubit, di);
            cache->completion_queue.pop_front();
        }
    }
}

void
GulliverSimulator::tick_prefetch()  {
    uint& ai = major_detector_register;
    uint& aj = minor_detector_table[ai];

    if (ai >= detector_vector_register.size()-1) {
        // We only transition state if all DRAM requests
        // have been processed.
        if (dram_await_array.size() == 0) {
            update_state();
        }
        return;
    }
    
    uint di = detector_vector_register[ai];
    uint dj = detector_vector_register[aj];

    if (di > curr_max_detector || dj > curr_max_detector) {
        return;
    }

    addr_t address = to_address(di, dj, base_address, n_detectors);
    
    access(address, false);  // We don't care if there is a hit or miss.
    // Update ai, aj.
    aj++;
    if (aj >= detector_vector_register.size()) {
        ai++;
        if (minor_detector_table.count(ai) == 0) {
            minor_detector_table[ai] = ai + 1;
        }
    }
}

void
GulliverSimulator::tick_bfu() {
    // We simulate the stages of the pipeline in reverse order.
    // There are bfu_hw_threshold stages. The first stage
    // is the fetch stage, the rest are computational stages.
    bfu_idle = true;
    for (uint stage = bfu_hw_threshold-1; stage > 0; stage--) {
        tick_bfu_compute(stage-1);  // Input is the compute stage number
                                    // so the second stage of the pipeline
                                    // would be stage 0 (as it is the first
                                    // compute stage.
    }
    tick_bfu_fetch();
    if (bfu_idle) {
        update_state();
    }
}

void
GulliverSimulator::tick_bfu_compute(uint stage) {
    // Get data from the pipeline latch.
    for (uint f = 0; f < bfu_fetch_width; f++) {
        BFUPipelineLatch& latch = bfu_pipeline_latches[f][stage];
        if (!latch.valid) {
            continue;
        }

        if (latch.stalled || latch.proposed_matches.size() == 0) {
            if (stage < bfu_hw_threshold-2) {
                // Stall next pipeline stage as well.
                bfu_pipeline_latches[f][stage+1].stalled = true;
            }
            continue;
        }
        bfu_idle = false;

        // Add entry from proposed matches to the running matching.
        std::map<uint, uint> matching(latch.base_entry.running_matching);
        fp_t matching_weight = latch.base_entry.matching_weight;
        auto match = latch.proposed_matches.top();
        // Update matching.
        uint di = detector_vector_register[latch.base_entry.next_unmatched_index];
        uint dj = match.first;
        matching[di] = dj;
        matching[dj] = di;
        // Update weight
        matching_weight += match.second;
        // Compute next unmatched detector.
        uint next_unmatched_index = detector_vector_register.size();
        for (uint ai = latch.base_entry.next_unmatched_index + 1; 
                ai < detector_vector_register.size(); ai++) 
        {
            if (matching.count(detector_vector_register[ai]) == 0) {
                next_unmatched_index = ai;
                break;
            } 
        }
        // Push this result onto the stack if there is a next unmatched
        // detector.
        StackEntry ei = {matching, matching_weight, next_unmatched_index};
        BFUPipelineLatch& next = bfu_pipeline_latches[f][stage+1];
        if (next_unmatched_index < detector_vector_register.size()) {
            hardware_stack.push(ei);
            // Update next latch.
            if (stage < bfu_hw_threshold-2) {
                next.stalled = false;
                next.valid = true;
                next.base_entry = latch.base_entry;
                next.proposed_matches = latch.proposed_matches;
                next.proposed_matches.pop();
            }
        } else {
            if (matching_weight < best_matching_register.matching_weight) { 
                best_matching_register = ei;
            }

            if (stage < bfu_hw_threshold-2) {
                next.stalled = true;
            }
        }
    }
}

void
GulliverSimulator::tick_bfu_fetch() {
    for (uint f = 0; f < bfu_fetch_width; f++) {
        // Determine what reads need to be made to the register file.
        // As we don't ever write to the register file outside of PREFETCH,
        // we can support many parallel reads.

        BFUPipelineLatch& next = bfu_pipeline_latches[f][0];
        if (hardware_stack.empty()) {
            next.stalled = true;
            continue;
        }
        bfu_idle = false;

        bool stall_next = false;

        std::stack<std::pair<uint, fp_t>> proposed_matches;
        StackEntry ei = hardware_stack.top();
        uint di = detector_vector_register[ei.next_unmatched_index];

        std::vector<std::pair<uint, fp_t>> match_list;
        fp_t min_weight = std::numeric_limits<fp_t>::max();
        for (uint aj = ei.next_unmatched_index+1; 
                aj < detector_vector_register.size(); aj++) 
        {
            uint dj = detector_vector_register[aj];
            if (ei.running_matching.count(dj)) {
                continue;
            }
            addr_t address = to_address(di, dj, base_address, n_detectors);
            bool is_hit = access(address, true);
            // If there is no hit on the register file, then we have to go
            // to DRAM. We will stall the pipeline until DRAM gives us
            // the value.
            auto di_dj = std::make_pair(di, dj);
            if (is_hit) {
                // Get weight and add data to proposed_matches.
                fp_t w = path_table[di_dj].distance;
                match_list.push_back(std::make_pair(dj, w));
                if (w < min_weight) {
                    min_weight = w;
                }
            } else {
                stall_next = true;
            } 
        }

        if (!stall_next) {
            for (auto p : match_list) {
                if (detector_vector_register.size() <= FILTER_CUTOFF 
                        || p.second < 1.2 * min_weight) 
                {
                    proposed_matches.push(p);
                }
            }
        }
        // Update the latch data.
        next.valid = true;
        next.proposed_matches = proposed_matches;
        next.base_entry = ei;
        next.stalled = stall_next;
        if (!stall_next) {
            hardware_stack.pop();
        }
    }
}

bool
GulliverSimulator::access(addr_t address, bool set_evictable_on_hit) {
    // First, check that the data isn't already in the register file.
    // or in the replacement queue.
    bool is_hit = false;
    for (Register& r : register_file) {
        if (r.valid && r.address == address) {
            is_hit = true;
            r.last_use = 0;
            break;
        }
    }
    for (addr_t x : replacement_queue) {
        if (x == address) {
            is_hit = true;
            break;
        }
    }
    // Check the cache.
    if (cache != nullptr) {
        auto di_dj = from_address(address, base_address, n_detectors);
        uint di = di_dj.first;
        uint dj = di_dj.second;
        QubitCache::LoadStatus s = cache->access(curr_qubit, di, dj);
        if (s == QubitCache::LoadStatus::hit) {
            replacement_queue.push_back(address);  // Add to register file.
            return true;
        } else if (s == QubitCache::LoadStatus::unknown) {
            // We exit this function immediately because we don't
            // want to go to DRAM yet. We are still waiting on the
            // result.
            return false;
        } // Otherwise, go perform the DRAM access.
    }

    if (!is_hit) {
        // If there is no hit, also check that the request has not
        // already been put out.
        bool is_waiting = false;
        for (auto x : dram_await_array) {
            if (address == x) {
                is_waiting = true;
                break;
            }
        }
        // Send out transaction to DRAM.
        if (!is_waiting && dram->WillAcceptTransaction(address, false)) {
            dram->AddTransaction(address, false);
            dram_await_array.push_back(address);
#ifdef GSIM_DEBUG
            std::cout << "\t[Gulliver] Requested " << std::hex
                << address << std::dec << "\n";
#endif
        }
    }

    return is_hit;
}

void
GulliverSimulator::update_state() {
    switch (state) {
    case State::ccomp:
#ifdef GSIM_DEBUG
        std::cout << "PREFETCH\n";
#endif
        state = State::prefetch;
        break;
    case State::prefetch:
#ifdef GSIM_DEBUG
        std::cout << "BFU\n";
#endif
        state = State::bfu;
        hardware_stack.push((StackEntry) {std::map<uint, uint>(), 0.0, 0});
        break;
    case State::bfu:
#ifdef GSIM_DEBUG
        std::cout << "IDLE\n";
#endif
        state = State::idle;
        break;
    default:
        break;
    }
}

void
GulliverSimulator::clear() {
    dram_await_array.clear();
    detector_vector_register.clear();
    
    ccomp_detector_register = 0;

    major_detector_register = 0;
    minor_detector_table.clear();
    minor_detector_table[0] = 1;

    while (!hardware_stack.empty()) {
        hardware_stack.pop();
    }
    for (uint f = 0; f < bfu_fetch_width; f++) {
        for (auto& latch : bfu_pipeline_latches[f]) {
            latch.stalled = false;
            latch.valid = false;
        }
    }
    best_matching_register = {
        std::map<uint, uint>(), 
        std::numeric_limits<fp_t>::max(),
        0 
    };
    replacement_queue.clear();
    if (cache == nullptr) {
        state = State::prefetch;
#ifdef GSIM_DEBUG
        std::cout << "PREFETCH\n";
#endif
    } else {
        state = State::ccomp;
#ifdef GSIM_DEBUG
        std::cout << "CCOMP\n";
#endif
    }
}

addr_t get_base_address(uint8_t bankgroup, uint8_t bank, uint32_t row_offset,
        dramsim3::Config * config) 
{
    addr_t base_address = bankgroup << config->bg_pos
                        | bank << config->ba_pos
                        | row_offset << config->ro_pos;
    base_address <<= config->shift_bits;
    return base_address;
}

addr_t
to_address(uint di, uint dj, addr_t base, uint n_detectors) {
    uint bdi = bound_detector(di, n_detectors);
    uint bdj = bound_detector(dj, n_detectors);
    addr_t x = base + (bdi*n_detectors + bdj) * sizeof(fp_t);
    return x;
}

std::pair<uint, uint>
from_address(addr_t x, addr_t base, uint n_detectors) {
    x -= base;
    uint bdi = x / (sizeof(fp_t) * n_detectors);
    uint bdj = x / sizeof(fp_t) - bdi*n_detectors;
    uint di = unbound_detector(bdi, n_detectors);
    uint dj = unbound_detector(bdj, n_detectors);
    std::pair<uint, uint> di_dj = std::make_pair(di, dj);
    return di_dj;
}

uint
bound_detector(uint d, uint n_detectors) {
    return d == BOUNDARY_INDEX ? n_detectors-1 : d;
}

uint
unbound_detector(uint d, uint n_detectors) {
    return d == n_detectors-1 ? BOUNDARY_INDEX : d;
}

}  // gulliver
}  // qrc
