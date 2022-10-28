/*
 *  author: Suhas Vittal
 *  date:   7 September 2022
 * */

#include "hyperion/simulator.h"
#include <limits>

namespace qrc {
namespace hyperion {

static const uint BFU_SORT_STAGES = 4;
static const uint RADIX_WIDTH = 16 / (BFU_SORT_STAGES);

HyperionSimulator::HyperionSimulator(dramsim3::MemorySystem * dram, 
        std::map<addr_t, bool> * memory_event_table,
        const PathTable& path_table,
        const HyperionSimulatorParams& params)
    :
    /* Statistics */
    prefetch_cycles(0),
    bfu_cycles(0),
    /* Microarchitecture */
    dram(dram),
    register_file(params.n_registers, (Register){0x0, 0, false}),
    dram_await_array(),
    detector_vector_register(),
    mean_weight_register(0),
    access_counter(0),
    major_detector_register(0),
    minor_detector_table(),
    hardware_deque(),
    bfu_pipeline_latches(params.bfu_fetch_width, 
            std::vector<BFUPipelineLatch>(
                BFU_SORT_STAGES+params.bfu_compute_stages)),
    best_matching_register(),
    replacement_queue(),
    state(HyperionSimulator::State::prefetch),
    bfu_idle(false),
    /* Data */
    memory_event_table(memory_event_table),
    path_table(path_table),
    /* Config parameters */
    curr_max_detector(0),
    n_detectors(params.n_detectors),
    n_detectors_per_round(params.n_detectors_per_round),
    bfu_fetch_width(params.bfu_fetch_width),
    bfu_compute_stages(params.bfu_compute_stages),
    curr_qubit(0),
    base_address(0),
    has_boundary(false)
{
    clear();
    load_base_address(params.bankgroup, params.bank, params.row_offset);
}

bool
HyperionSimulator::load_detectors(const std::vector<uint>& detector_array) {
    if (detector_array.size() == 0) {
        state = State::idle;
        return false;
    }
    clear();
    detector_vector_register = detector_array;
#ifdef HSIM_DEBUG
    std::cout << "LOAD:";
    for (uint d : detector_array) {
        std::cout << " " << d; 
    }
    std::cout << "\n";
#endif
    has_boundary = detector_vector_register.back() == BOUNDARY_INDEX;
    // Compute the critical index.
    curr_max_detector = n_detectors_per_round;

    return true;
}

void
HyperionSimulator::load_path_table(const PathTable& pt) {
    path_table = pt;  
}

void
HyperionSimulator::load_qubit_number(uint qubit) {
    curr_qubit = qubit;
}

void
HyperionSimulator::load_base_address(uint8_t bg, uint8_t ba, uint32_t ro_offset) {
    if (dram == nullptr) {
        base_address = 0x0;
    } else {
        base_address = get_base_address(bg, ba, ro_offset, dram->config_);
    }
}

void
HyperionSimulator::tick() {
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
#ifdef HSIM_DEBUG
        std::cout << "\t[Hyperion] Installed " << std::hex 
            << address << std::dec << "\n";
#endif
    }
    // Check if any transactions to DRAM have completed.
    if (dram != nullptr) {
        for (auto it = dram_await_array.begin(); it != dram_await_array.end(); ) {
            addr_t address = *it;
            if (memory_event_table->count(address) 
                    && memory_event_table->at(address)) 
            {
                it = dram_await_array.erase(it);
                memory_event_table->at(address) = false;
                // Add to replacement queue.
                auto di_dj = from_address(address, base_address, n_detectors);
#ifdef HSIM_DEBUG
                std::cout << "\t[DRAM] Retrieved " << std::hex
                    << address << std::dec << "(" << di_dj.first
                    << "," << di_dj.second << ")\n";
#endif
                uint di = di_dj.first;
                replacement_queue.push_back(address);
            } else {
                it++;
            }
        }
    }

    switch (state) {
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
HyperionSimulator::sig_end_round(uint rounds_ended) {
    // We subtract by 1 because we don't count the boundary.
    curr_max_detector += n_detectors_per_round * rounds_ended;
    if (curr_max_detector > n_detectors-1) {
        curr_max_detector = n_detectors - 1;
    }
    major_detector_register = 0;
#ifdef HSIM_DEBUG
    std::cout << "SIGNAL -- ROUND END | Max detector = " 
        << curr_max_detector << "\n";
#endif
}

bool
HyperionSimulator::is_idle() {
    return state == State::idle;
}

void
HyperionSimulator::force_idle() {
    state = State::idle;
}

std::map<uint, uint>
HyperionSimulator::get_matching() {
    return best_matching_register.running_matching;
}

void
HyperionSimulator::reset_stats() {
    prefetch_cycles = 0;
    bfu_cycles = 0;
}

uint64_t
HyperionSimulator::rowhammer_flips() {
    if (dram == nullptr) return 0;
    return dram->dram_system_->rowhammer_flips();
}

uint64_t
HyperionSimulator::row_activations() {
    if (dram == nullptr) return 0;
    return dram->dram_system_->row_activations();
}

void
HyperionSimulator::tick_prefetch()  {
    uint& ai = major_detector_register;
    if (minor_detector_table.count(ai) == 0) {
        // Start with boundary first if it exists.
        if (has_boundary) {
            minor_detector_table[ai] = detector_vector_register.size() - 1;
        } else {
            minor_detector_table[ai] = ai + 1;
        }
    }
    uint& aj = minor_detector_table[ai];
    
    while (aj >= detector_vector_register.size()
        && ai < detector_vector_register.size()-1)
    {
        ai++;
        aj = minor_detector_table[ai];
    }

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

//    if (curr_max_detector > WINDOW_SIZE*n_detectors_per_round
//            && detector_vector_register.size() > FILTER_CUTOFF
//            && di < (curr_max_detector - WINDOW_SIZE*n_detectors_per_round)) 
//    {
//        ai++;
//        return;
//    }

    if (di != BOUNDARY_INDEX && di > curr_max_detector) {
        return;
    }

    if (dj != BOUNDARY_INDEX && dj > curr_max_detector) {
        uint next_di = detector_vector_register[ai+1];
        if (next_di == BOUNDARY_INDEX || next_di < curr_max_detector) {
            ai++;
        }
        return;
    }

    addr_t address = to_address(di, dj, base_address, n_detectors);
    
    access(address, false);  // We don't care if there is a hit or miss.
    fp_t raw_weight = path_table[std::make_pair(di, dj)].distance;
    uint32_t w = (uint32_t) (raw_weight * MWPM_INTEGER_SCALE);
    mean_weight_register += w;
    access_counter++;
    // Update ai, aj.
    // Update varies on whether or not boundary is in the DVR.
    if (has_boundary && dj == BOUNDARY_INDEX) {
        aj = ai + 1;  // Reset, we have finished the boundary.
    } else {
        aj++;
    }

    if ((has_boundary && aj >= detector_vector_register.size()-1) 
            || (!has_boundary && aj >= detector_vector_register.size())) 
    {
        aj = detector_vector_register.size();
        ai++;
    }
}

void
HyperionSimulator::tick_bfu() {
    // We simulate the stages of the pipeline in reverse order.
    // There is one FETCH stage, four SORT stages, and 
    // a variable number of COMPUTE stages.
    bfu_idle = true;
    for (uint stage = bfu_compute_stages; stage > 0; stage--) {
        tick_bfu_compute(stage-1);  // Input is the compute stage number
                                    // so the second stage of the pipeline
                                    // would be stage 0 (as it is the first
                                    // compute stage.
    }
    for (uint stage = BFU_SORT_STAGES; stage > 0; stage--) {
        tick_bfu_sort(stage-1);
    }
    tick_bfu_fetch();
    if (bfu_idle) {
        update_state();
    }
}

void
HyperionSimulator::tick_bfu_compute(uint stage) {
    const uint STAGE_OFFSET = BFU_SORT_STAGES;
    // Get data from the pipeline latch.
    for (uint f = 0; f < bfu_fetch_width; f++) {
        BFUPipelineLatch& latch = bfu_pipeline_latches[f][stage+STAGE_OFFSET];
        if (!latch.valid) {
            continue;
        }

        if (latch.stalled || latch.proposed_matches.size() == 0) {
            if (stage < bfu_compute_stages-1) {
                // Stall next pipeline stage as well.
                bfu_pipeline_latches[f][stage+STAGE_OFFSET+1].stalled = true;
            }
            continue;
        }
        bfu_idle = false;

        // Add entry from proposed matches to the running matching.
        std::map<uint, uint> matching(latch.base_entry.running_matching);
        uint32_t matching_weight = latch.base_entry.matching_weight;
        auto match = latch.proposed_matches.front();
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
        DequeEntry ei = {matching, matching_weight, next_unmatched_index};
        BFUPipelineLatch& next = bfu_pipeline_latches[f][stage+STAGE_OFFSET+1];
        if (next_unmatched_index < detector_vector_register.size()) {
#ifdef MATCHING_PQ
            hardware_deque.push(ei);
#else
            hardware_deque.push_back(ei);
#endif
            // Update next latch.
            if (stage < bfu_compute_stages-1) {
                next.stalled = false;
                next.valid = true;
                next.base_entry = latch.base_entry;
                next.proposed_matches = latch.proposed_matches;
                next.proposed_matches.pop_front();
            }
        } else {
            if (matching_weight < best_matching_register.matching_weight) { 
                best_matching_register = ei;
            }

            if (stage < bfu_compute_stages-1) {
                next.stalled = true;
            }
        }
    }
}

void
HyperionSimulator::tick_bfu_sort(uint stage) {
    // Get pipeline latch. 
    for (uint f = 0; f < bfu_fetch_width; f++) {
        BFUPipelineLatch& latch = bfu_pipeline_latches[f][stage];
        if (!latch.valid) {
            continue;
        }

        if (latch.stalled) {
            bfu_pipeline_latches[f][stage+1].stalled = true; 
            continue;
        }
        bfu_idle = false;
        // Get radix bins from prior round.
        auto proposed_matches = latch.proposed_matches;
        uint8_t nibble_pos = stage;
        std::array<std::vector<uint>, 1<<RADIX_WIDTH> bins;
        for (uint i = 0; i < proposed_matches.size(); i++) {
            uint32_t w = proposed_matches[i].second;
            // Get nibble (4-bits) at relevant position.
            uint8_t nibble = (w >> (nibble_pos << 2)) & 0xf;
            bins[nibble].push_back(i);
        }
        // Merge into one array.
        std::deque<std::pair<uint, uint32_t>> new_matches;
        for (auto bucket : bins) {
            for (auto i : bucket) {
                new_matches.push_back(proposed_matches[i]); 
            }
        }

        BFUPipelineLatch& next = bfu_pipeline_latches[f][stage+1];
        next.proposed_matches = new_matches;
        next.base_entry = latch.base_entry;
        next.stalled = false;
        next.valid = true;
    }
}

void
HyperionSimulator::tick_bfu_fetch() {
    for (uint f = 0; f < bfu_fetch_width; f++) {
        // Determine what reads need to be made to the register file.
        // As we don't ever write to the register file outside of PREFETCH,
        // we can support many parallel reads.

        BFUPipelineLatch& next = bfu_pipeline_latches[f][0];
        if (hardware_deque.empty()) {
            next.stalled = true;
            continue;
        }
        bfu_idle = false;

        bool stall_next = false;

        std::deque<std::pair<uint, uint32_t>> proposed_matches;
#ifdef MATCHING_STACK
        DequeEntry ei = hardware_deque.back();
#elifdef MATCHING_FIFO
        DequeEntry ei = hardware_deque.front();
#else
        DequeEntry ei = hardware_deque.top();
#endif
        uint di = detector_vector_register[ei.next_unmatched_index];

        std::vector<std::pair<uint, uint32_t>> match_list;
        uint32_t min_weight = std::numeric_limits<uint32_t>::max();
        for (uint aj = ei.next_unmatched_index+1; 
                aj < detector_vector_register.size(); aj++) 
        {
            uint dj = detector_vector_register[aj];
            if (ei.running_matching.count(dj)) {
                continue;
            }
//            if (dj != BOUNDARY_INDEX && di != BOUNDARY_INDEX
//                    && detector_vector_register.size() > FILTER_CUTOFF
//                    && (dj-di) > WINDOW_SIZE*n_detectors_per_round)
//            {
//                continue;
//            }
            addr_t address = to_address(di, dj, base_address, n_detectors);
            bool is_hit = access(address, true);
            // If there is no hit on the register file, then we have to go
            // to DRAM. We will stall the pipeline until DRAM gives us
            // the value.
            auto di_dj = std::make_pair(di, dj);
            if (is_hit) {
                // Get weight and add data to proposed_matches.
                fp_t raw_weight = path_table[di_dj].distance;
                // Convert to integer (should be like this in DRAM).
                uint32_t w = (uint32_t)(raw_weight * MWPM_INTEGER_SCALE);
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
                uint dj = p.first;
                uint32_t w = p.second;
                if (detector_vector_register.size() <= FILTER_CUTOFF
                    || w <= mean_weight_register
                    || w == min_weight)
                {
                    proposed_matches.push_back(p);
                }
            }
        }
        // Update the latch data.
        next.valid = true;
        next.proposed_matches = proposed_matches;
        next.base_entry = ei;
        next.stalled = stall_next;
        if (!stall_next) {
#ifdef MATCHING_STACK
            hardware_deque.pop_back();
#elifdef MATCHING_FIFO
            hardware_deque.pop_front();
#else
            hardware_deque.pop();
#endif
        }
    }
}

bool
HyperionSimulator::access(addr_t address, bool set_evictable_on_hit) {
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

    if (!is_hit) {
        // If there is no hit, also check that the request has not
        // already been put out.
        if (dram != nullptr) {
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
#ifdef HSIM_DEBUG
                std::cout << "\t[Hyperion] Requested " << std::hex
                    << address << std::dec << "\n";
#endif
            }
        } else {
            // Just directly pull the data from SRAM and add it to the
            // replacement queue.
            replacement_queue.push_back(address);
        }
    }

    return is_hit;
}

void
HyperionSimulator::update_state() {
    switch (state) {
    case State::prefetch:
        mean_weight_register /= access_counter;
#ifdef HSIM_DEBUG
        std::cout << "[Mean Weight] " << mean_weight_register << "\n";
        std::cout << "BFU\n";
#endif
        state = State::bfu;
#ifdef MATCHING_PQ
        hardware_deque.push((DequeEntry) {std::map<uint, uint>(), 0, 0});
#else
        hardware_deque.push_back((DequeEntry) {std::map<uint, uint>(), 0, 0});
#endif
        break;
    case State::bfu:
#ifdef HSIM_DEBUG
        std::cout << "IDLE\n";
        std::cout << "\tcycles in prefetch: " << prefetch_cycles << "\n";
        std::cout << "\tcycles in bfu: " << bfu_cycles << "\n";
        std::cout << "\thamming weight: " << detector_vector_register.size() << "\n";
#endif
        state = State::idle;
        break;
    default:
        break;
    }
}

void
HyperionSimulator::clear() {
    dram_await_array.clear();
    detector_vector_register.clear();
    mean_weight_register = 0;
    access_counter = 0;
    
    major_detector_register = 0;
    minor_detector_table.clear();

#ifdef MATCHING_PQ
    while (!hardware_deque.empty()) {
        hardware_deque.pop();
    }
#else
    hardware_deque.clear();
#endif
    for (uint f = 0; f < bfu_fetch_width; f++) {
        for (auto& latch : bfu_pipeline_latches[f]) {
            latch.stalled = false;
            latch.valid = false;
        }
    }
    best_matching_register = {
        std::map<uint, uint>(), 
        std::numeric_limits<uint32_t>::max(),
        0 
    };
    replacement_queue.clear();
    state = State::prefetch;
#ifdef HSIM_DEBUG
    std::cout << "(clear)PREFETCH\n";
#endif
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
    addr_t x = base + (bdi*n_detectors + (bdj-bdi-1)) * sizeof(uint32_t);
    return x;
}

std::pair<uint, uint>
from_address(addr_t x, addr_t base, uint n_detectors) {
    x -= base;
    uint bdi = x / (sizeof(uint32_t) * n_detectors);
    uint bdj = x / sizeof(uint32_t) - bdi*n_detectors + bdi + 1;
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

}  // hyperion
}  // qrc
