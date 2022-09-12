/*
 *  author: Suhas Vittal
 *  date:   7 September 2022
 * */

#include "gulliver/simulator.h"
#include <limits>

//#define GSIM_DEBUG

GulliverSimulator::GulliverSimulator(dramsim3::MemorySystem * dram, 
        std::map<addr_t, bool> * memory_event_table,
        std::map<uint, std::vector<uint>> * adjacency_lists,
        const PathTable& path_table,
        const GulliverSimulatorParams& params)
    :
    /* Statistics */
    prefetch_cycles(0),
    bfu_cycles(0),
    /* Microarchitecture */
    dram(dram),
    register_file(params.n_registers),
    next_dram_request_register(std::make_pair(0,0)),
    dram_await_array(),
    detector_vector_register(),
    predecode_scoreboard(params.bfu_hw_threshold),
    index_register(0),
    neighbor_register(0),
    hardware_stack(),
    bfu_pipeline_latches(params.bfu_fetch_width, 
            std::vector<BFUPipelineLatch>(params.bfu_hw_threshold-1)),
    best_matching_register(),
    replacement_queue(),
    state(GulliverSimulator::State::prefetch),
    bfu_idle(false),
    /* Config parameters */
    memory_event_table(memory_event_table),
    path_table(path_table),
    adjacency_lists(adjacency_lists),
    n_detectors(params.n_detectors),
    bfu_fetch_width(params.bfu_fetch_width),
    bfu_hw_threshold(params.bfu_hw_threshold),
    base_address(0)
{
    load_base_address(params.bankgroup, params.bank, params.row_offset);
}

bool
GulliverSimulator::load_detectors(const std::vector<uint>& detector_array) {
    clear();
    detector_vector_register = detector_array;
#ifdef GSIM_DEBUG
    std::cout << "PREDECODE\n";
#endif
    return detector_array.size() <= bfu_hw_threshold;
}

void
GulliverSimulator::load_path_table(const PathTable& pt) {
    clear();
    path_table = pt;  
}

void
GulliverSimulator::load_base_address(uint8_t bg, uint8_t ba, uint32_t ro_offset) {
    clear();
    base_address = get_base_address(bg, ba, ro_offset, dram->config_);
}

void
GulliverSimulator::tick() {
    dram->ClockTick();
    
    // Go perform any replacements if necessary.
    if (replacement_queue.size()) {
        addr_t address = replacement_queue.front();
        // Check if there are any evictable registers.
        bool is_installed = false;
        for (Register& r : register_file) {
            if (!r.valid || r.evictable) {
                r.address = address; 
                r.evictable = false;
                r.valid = true;
                is_installed = true;
                break;
            }
        }
        
        if (is_installed) {
            replacement_queue.pop_front();
#ifdef GSIM_DEBUG
            std::cout << "\t[Gulliver] Installed " << std::hex 
                << address << std::dec << "\n";
#endif
        }
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
#ifdef GSIM_DEBUG
            std::cout << "\t[DRAM] Retrieved " << std::hex
                << address << std::dec << "\n";
#endif
            replacement_queue.push_back(address);
        } else {
            it++;
        }
    }

    switch (state) {
    case GulliverSimulator::State::prefetch:
        prefetch_cycles++;
        tick_prefetch();
        break;
    case GulliverSimulator::State::predecode:
        tick_predecode();
        break;
    case GulliverSimulator::State::consolidate:
        tick_consolidate();
        break;
    case GulliverSimulator::State::bfu:
        bfu_cycles++;
        tick_bfu();
        break;
    default:
        break;
    }
}

bool
GulliverSimulator::is_idle() {
    return state == GulliverSimulator::State::idle;
}

std::map<uint, uint>
GulliverSimulator::get_matching() {
    return best_matching_register.running_matching;
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
GulliverSimulator::tick_predecode() {
    if (index_register < detector_vector_register.size()) {
        uint di = detector_vector_register[index_register];
        std::vector<uint> neighbors = adjacency_lists->at(di);
        // Check if neighbor is in in the 
        // DVR (detector vector register).
        bool in_dvr = false;
        uint dj = neighbors[neighbor_register];
        for (uint d : detector_vector_register) {
            if (dj == d) {
                in_dvr = true;
                break;
            }
        }
        // If in the DVR, then update entry in the scoreboard.
        // Also, if di >= dj, then we have already seen this combination,
        // so we can just skip if that's the case.
        if (index_register < neighbor_register && in_dvr) {
            addr_t weight_addr = to_address(di, dj, base_address, n_detectors);
            fp_t w = path_table[std::make_pair(di, dj)].distance;
            bool is_hit = access(weight_addr, false);
            if (!is_hit) {
                return;
            }

            PDScoreboardEntry& ei = predecode_scoreboard[index_register];
            PDScoreboardEntry& ej = predecode_scoreboard[neighbor_register];
            if (ei.mate_weight > w) {
                ei.best_mate_index = neighbor_register;
                ej.n_suitors++;
                ei.mate_weight = w;
            }
            if (ej.mate_weight > w) {
                ej.best_mate_index = index_register;
                ei.n_suitors++;
                ej.mate_weight = w;
            }
        }
        // Update index and neighbor registers.
        neighbor_register++;
        if (neighbor_register >= neighbors.size()) {
            index_register++;
            neighbor_register = 0;
        }
    } else {
        state = GulliverSimulator::State::consolidate;
        index_register = 0;
#ifdef GSIM_DEBUG
        std::cout << "CONSOLIDATE\n";
#endif
    }
}

void
GulliverSimulator::tick_consolidate() {
    // We use the best matching register to consolidate results
    // from predecode.
    //
    // As we exit the consolidate state and enter prefetch,
    // we will reset the best matching register and push the contents
    // onto the hardware stack used in the BFU pipeline.
    std::map<uint, uint>& init_matching = best_matching_register.running_matching;
    fp_t& init_weight = best_matching_register.matching_weight;
    if (index_register < detector_vector_register.size()) {
        // Retrieve data
        uint di = detector_vector_register[index_register];
        if (init_matching.count(di)) {
            index_register++;   // Auto-increment if we have already assigned di.
        }
        PDScoreboardEntry& ei = predecode_scoreboard[index_register];
        // Retrieve neighbor data. 
        uint j = ei.best_mate_index;
        uint dj = detector_vector_register[j]; 
        PDScoreboardEntry& ej = predecode_scoreboard[j];
        // Make sure dj "likes" di.
        if (index_register == ej.best_mate_index
                && ei.n_suitors == 1 && ej.n_suitors == 1)
        {
            init_matching[di] = dj;
            init_matching[dj] = di;
            init_weight += ei.mate_weight;
        }
        index_register++;
    } else {
        StackEntry init = {
            init_matching,
            init_weight,
            0
        };
        // Remove all entries in init_matching from
        // the DVR.
        for (auto it = detector_vector_register.begin();
                it != detector_vector_register.end(); )
        {
            if (init_matching.count(*it)) {
                it = detector_vector_register.erase(it);
            } else {
                it++;
            }
        }
        if (detector_vector_register.empty()) {
            state = GulliverSimulator::State::idle;
        } else {
            // Reset best_matching_register
            best_matching_register.running_matching = std::map<uint, uint>();
            best_matching_register.matching_weight = 0.0;
            hardware_stack.push(init);

            state = GulliverSimulator::State::prefetch;
#ifdef GSIM_DEBUG
            std::cout << "PREFETCH\n";
#endif
        }
    }
}

void
GulliverSimulator::tick_prefetch()  {
    uint& ai = next_dram_request_register.first;
    uint& aj = next_dram_request_register.second;

    if (ai >= detector_vector_register.size()-1) {
        // We only transition state if all DRAM requests
        // have been processed.
        if (dram_await_array.size() == 0) {
            state = GulliverSimulator::State::bfu;
#ifdef GSIM_DEBUG
            std::cout << "BFU\n";
#endif
        }
        return;
    }
    
    uint di = detector_vector_register[ai];
    uint dj = detector_vector_register[aj];
    addr_t address = to_address(di, dj, base_address, n_detectors);
    std::pair<uint, uint> di_dj = std::make_pair(di, dj);
    
    access(address, false);  // We don't care if there is a hit or miss.
    // Update di, dj.
    aj++;
    if (aj >= detector_vector_register.size()) {
        ai++;
        aj = ai+1;
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
        state = GulliverSimulator::State::idle;
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
                proposed_matches.push(std::make_pair(dj, w));
            } else {
                stall_next = true;
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
            r.evictable = set_evictable_on_hit;
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
GulliverSimulator::clear() {
    for (Register& r : register_file) {
        r.evictable = true;
    }    
    next_dram_request_register = std::make_pair(0,1);
    dram_await_array.clear();
    detector_vector_register.clear();
    for (PDScoreboardEntry& ei : predecode_scoreboard) {
        ei.best_mate_index = 0;
        ei.n_suitors = 0;
        ei.mate_weight = std::numeric_limits<fp_t>::max();
    }
    index_register = 0;
    neighbor_register = 0;
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
    state = GulliverSimulator::State::predecode;
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
    addr_t x = base + (bdi*n_detectors + bdj-bdi-1) * sizeof(fp_t);
    return x;
}

std::pair<uint, uint>
from_address(addr_t x, addr_t base, uint n_detectors) {
    x -= base;
    uint bdi = x / (sizeof(fp_t) * n_detectors);
    uint bdj = x / sizeof(fp_t) - bdi*n_detectors + bdi + 1;
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
