/*
 *  author: Suhas Vittal
 *  date:   7 September 2022
 * */

#include "gulliver/simulator.h"

GulliverSimulator::GulliverSimulator(dramsim3::MemorySystem * dram, 
        std::map<std::pair<uint, uint>, bool> * memory_event_table,
        const std::map<std::pair<uint, uint>, fp_t>& weight_table,
        const GulliverSimulatorParams& params)
    :dram(dram),
    register_file(params.n_registers),
    next_dram_request_register(std::make_pair(0,0)),
    dram_await_array(),
    detector_vector_register(),
    hardware_stack(),
    pipeline_latches(params.bfu_fetch_width, 
            std::vector<PipelineLatch>(params.bfu_hw_threshold-1)),
    best_matching_register(),
    replacement_queue(),
    state(GulliverSimulator::State::prefetch),
    bfu_idle(false),
    memory_event_table(memory_event_table),
    weight_table(weight_table),
    n_detectors(params.n_detectors),
    bfu_fetch_width(params.bfu_fetch_width),
    bfu_hw_threshold(params.bfu_hw_threshold),
    bankgroup(params.bankgroup),
    bank(params.bank),
    row_offset(params.row_offset),
    base_address(0)
{
    clear();
    // Setup base address. 
    base_address = get_base_address(bankgroup, bank, row_offset, dram->config_);
}

void
GulliverSimulator::load(const std::vector<uint>& detector_array) {
    clear();
    detector_vector_register = detector_array;
}

void
GulliverSimulator::tick() {
    dram->ClockTick();
    
    // Go perform any replacements if necessary.
    while (replacement_queue.size()) {
        addr_t address = replacement_queue.front();
        // Check if there are any evictable registers.
        bool is_installed = false;
        for (Register& r : register_file) {
            if (!r.valid || r.evictable) {
                r.address = address; 
                r.evictable = false;
                r.valid = true;
                is_installed = true;
            }
        }
        
        if (is_installed) {
            replacement_queue.pop_front();
        } else {
            break;
        }
    }
    // Check if any transactions to DRAM have completed.
    for (auto it = dram_await_array.begin(); it != dram_await_array.end(); ) {
        auto di_dj = *it;
        if (memory_event_table->count(di_dj) 
                && memory_event_table->at(di_dj)) 
        {
            it = dram_await_array.erase(it);
            memory_event_table->at(di_dj) = false;
            // Add to replacement queue.
            uint di = di_dj.first;
            uint dj = di_dj.second;
            addr_t address = 
                to_address(di, dj, base_address, n_detectors);
            replacement_queue.push_back(address);
        } else {
            it++;
        }
    }

    switch (state) {
    case GulliverSimulator::State::prefetch:
        tick_prefetch();
        break;
    case GulliverSimulator::State::predecode:
        tick_predecode();
        break;
    case GulliverSimulator::State::bfu:
        tick_bfu();
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

void
GulliverSimulator::tick_prefetch()  {
    uint& ai = next_dram_request_register.first;
    uint& aj = next_dram_request_register.second;

    if (ai >= detector_vector_register.size()-1) {
        // We only transition state if all DRAM requests
        // have been processed.
        if (dram_await_array.size() == 0) {
            state = GulliverSimulator::State::predecode;
        }
        return;
    }
    
    uint di = detector_vector_register[ai];
    uint dj = detector_vector_register[aj];
    addr_t address = to_address(di, dj, base_address, n_detectors);
    std::pair<uint, uint> di_dj = std::make_pair(di, dj);
    // Send out transaction to DRAM.
    if (dram->WillAcceptTransaction(address, false)) {
        dram->AddTransaction(address, false);
        dram_await_array.push_back(di_dj);
        // Update di, dj.
        aj++;
        if (aj >= detector_vector_register.size()) {
            ai++;
            aj = ai+1;
        }
    }
}

void
GulliverSimulator::tick_predecode() {
    state = GulliverSimulator::State::bfu;
    StackEntry init = {
        std::map<uint, uint>(),
        0.0,
        0
    };
    hardware_stack.push(init);
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
        PipelineLatch& latch = pipeline_latches[f][stage];
        if (!latch.valid) {
            continue;
        }

        if (latch.stalled || latch.proposed_matches.size() == 0) {
            if (stage < bfu_hw_threshold-2) {
                // Stall next pipeline stage as well.
                pipeline_latches[f][stage+1].stalled = true;
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
        StackEntry e = {matching, matching_weight, next_unmatched_index};
        PipelineLatch& next = pipeline_latches[f][stage+1];
        if (next_unmatched_index < detector_vector_register.size()) {
            hardware_stack.push(e);
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
                best_matching_register = e;
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
        PipelineLatch& next = pipeline_latches[f][0];
        if (hardware_stack.empty()) {
            next.stalled = true;
            continue;
        }
        bfu_idle = false;

        std::stack<std::pair<uint, fp_t>> proposed_matches;
        StackEntry e = hardware_stack.top();
        uint di = detector_vector_register[e.next_unmatched_index];
        for (uint aj = e.next_unmatched_index+1; 
                aj < detector_vector_register.size(); aj++) 
        {
            uint dj = detector_vector_register[aj];
            if (e.running_matching.count(dj)) {
                continue;
            }
            addr_t address = to_address(di, dj, base_address, n_detectors);
            bool is_hit = false;

            for (Register& r : register_file) {
                if (r.valid && r.address == address) {
                    is_hit = true;
                    break;
                }
            }
            // If there is no hit on the register file, then we have to go
            // to DRAM. We will stall the pipeline until DRAM gives us
            // the value.
            auto di_dj = std::make_pair(di, dj);
//          if (is_hit) {
                // Get weight and add data to proposed_matches.
                fp_t w = weight_table[di_dj];
                proposed_matches.push(std::make_pair(dj, w));
                next.stalled = false;
/*          } else {
                std::cout << 
                    "[GulliverSimulator] error: data not in register file for: " 
                    << di << "," << dj << ".\n";
                // Check if di_dj are in the dram_await_array
                // If they are not, then we send a transaction.
                bool is_waiting = false;
                for (auto ei_ej : dram_await_array) {
                    if (di_dj == ei_ej) {
                        is_waiting = true;
                        break;
                    }
                }
                if (!is_waiting && dram->WillAcceptTransaction(address, false)) {
                    dram->AddTransaction(address, false);
                    dram_await_array.push_back(di_dj);
                }
                latch.stalled = true;
            } 
*/      }
        // Update the latch data.
        next.valid = true;
        next.proposed_matches = proposed_matches;
        next.base_entry = e;
        if (!next.stalled) {
            hardware_stack.pop();
        }
    }
}

void
GulliverSimulator::clear() {
    for (Register& r : register_file) {
        r.address = 0x0;
        r.evictable = true;
        r.valid = false;
    }    
    next_dram_request_register = std::make_pair(0,1);
    dram_await_array.clear();
    detector_vector_register.clear();
    while (!hardware_stack.empty()) {
        hardware_stack.pop();
    }
    for (uint f = 0; f < bfu_fetch_width; f++) {
        for (PipelineLatch& latch : pipeline_latches[f]) {
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
    state = GulliverSimulator::State::prefetch;
}

addr_t get_base_address(uint8_t bankgroup, uint8_t bank, uint32_t row_offset,
        dramsim3::Config * config) 
{
    addr_t base_address = bankgroup << config->bg_pos
                        || bank << config->ba_pos
                        || row_offset << config->ro_pos;
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
