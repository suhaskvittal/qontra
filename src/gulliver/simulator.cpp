/*
 *  author: Suhas Vittal
 *  date:   27 October 2022
 * */

#include "gulliver/simulator.h"

namespace qrc {
namespace gulliver {

GulliverSimulator::GulliverSimulator(
        dramsim3::MemorySystem * dram,
        std::map<addr_t, bool> * memory_event_table,
        const DecodingGraph& decoding_graph,
        MWPMDecoder * decoder,
        uint search_width)
    :dram(dram),
    traverse_cycles(0),
    combine_cycles(0),
    max_active_syndromes_in_box(0),
    syndrome_register(),
    hamming_weight_register(0),
    current_qubit(0),
    traverse_pipeline_latches(1, 
                (TraversePipelineLatch) {
                std::vector<uint>(), std::map<uint, uint>(),
                0, false, false}
            ),
    cumulative_matching_register(),
    matchings(),
    traverse_idle(false),
    decoder(decoder),
    decoding_graph(decoding_graph),
    memory_event_table(memory_event_table),
    search_width(search_width)
{}

void
GulliverSimulator::load_syndrome(const std::vector<uint8_t>& syndrome) {
    clear();
    syndrome_register = syndrome;
    // Get hamming weight.
    hamming_weight_register = std::accumulate(
                                syndrome_register.begin(),
                                syndrome_register.end(),
                                0);
#ifdef GSIM_DEBUG
    std::cout << "LOAD: syndrome of weight " << hamming_weight_register
        << ", detectors =";
    for (uint i = 0; i < syndrome.size(); i++) {
        if (syndrome[i]) std::cout << " " << i;
    }
    std::cout << "\n";
#endif
    state = State::traverse;
}

void
GulliverSimulator::tick() {
    // DRAM should be ticked externally.
    // Check memory event table.
    for (auto it = dram_await_array.begin(); it != dram_await_array.end(); ) {
        addr_t x = *it;
        if (memory_event_table->count(x) && memory_event_table->at(x)) {
            // Perform MWPM on the corresponding syndrome.
            auto syndrome = addr_to_syndrome[x];
#ifdef GSIM_DEBUG
            std::cout << "\t[Gulliver] Received matching at address "
                        << std::hex << x << std::dec << ":";
#endif
            DecoderShotResult res = decoder->decode_error(syndrome);
            auto matching = res.matching;
            for (auto pair : matching) {
                uint di = pair.first;
                uint dj = pair.second;
#ifdef GSIM_DEBUG
                std::cout << " " << di << " --> " << dj;
#endif
                if (!matchings.count(di)) {
                    matchings[di] = std::set<uint>();
                }
                matchings[di].insert(dj);
            }
#ifdef GSIM_DEBUG
            std::cout << "\n";
#endif
            // Delete entry from await array and clear table entry.
            it = dram_await_array.erase(it);
            memory_event_table->at(x) = false;
        } else {
            it++;
        }
    }
    // Now, tick the rest of the microarchitecture.
    switch (state) {
    case State::traverse:
        tick_traverse();
        traverse_cycles++;
        break;
    case State::combine:
        tick_combine();
        combine_cycles++;
        break;
    case State::idle:
        break;
    }
}

bool
GulliverSimulator::is_idle() {
    return state == State::idle;
}

void
GulliverSimulator::force_idle() {
    state = State::idle;
}

std::map<uint, uint>
GulliverSimulator::get_matching() {
    return cumulative_matching_register;
}

void
GulliverSimulator::reset_stats() {
    traverse_cycles = 0;
    combine_cycles = 0;
    max_active_syndromes_in_box = 0;
}

void
GulliverSimulator::tick_traverse() {
    traverse_idle = true;
    tick_traverse_mem();
    tick_traverse_read();
    if (traverse_idle) {
        if (dram_await_array.empty()) {
            update_state();
        }
    }
}

void
GulliverSimulator::tick_combine() {
    std::deque<std::map<uint, uint>> matching_list; 
    // Multipath execution style -- get all valid
    // matchings (i.e. corrections whose errors correspond
    // generate the syndrome).
    //
    // Note that in actuality, we would be checking
    // which of these matchings correspond to corrections
    // that generate the original syndrome. But this is dumber
    // so oh well.
    //
    // TODO make this more fine grained.
    matching_list.push_back(std::map<uint, uint>());
    for (auto pair : matchings) {
        uint di = pair.first;
        auto dj_set = pair.second;

        std::deque<std::map<uint, uint>> next_list;
        while (!matching_list.empty()) {
            auto m = matching_list.front();
            matching_list.pop_front();
            if (m.count(di)) {
                next_list.push_back(m);
            } else {
#ifdef GSIM_DEBUG
                std::cout << "\tChecking matches for " << di << ":";
#endif
                for (uint dj : dj_set) {
#ifdef GSIM_DEBUG
                    std::cout << " " << dj;
#endif
                    if (m.count(dj)) {
                        next_list.push_back(m);
                        continue;
                    }
                    std::map<uint, uint> nm(m);
                    nm[di] = dj;
                    nm[dj] = di;
                    next_list.push_back(nm);
                }
#ifdef GSIM_DEBUG
                std::cout << "\n";
#endif
            }
        }
        matching_list = next_list;
    }
#ifdef GSIM_DEBUG
    std::cout << "\tNumber of matchings: " << matching_list.size() << "\n";
#endif

    // Return the first good matching.
#ifdef GSIM_DEBUG
    std::cout << "\tMatching sizes:";
#endif
    for (auto m : matching_list) {
#ifdef GSIM_DEBUG
        std::cout << " " << m.size();
#endif
        if (m.size() >= hamming_weight_register) {
            cumulative_matching_register = m;
            break;
        }
    }
#ifdef GSIM_DEBUG
    std::cout << "\n";
#endif

    update_state();
}

void
GulliverSimulator::tick_traverse_mem() {
    // Retrieve data from latch.
    TraversePipelineLatch& latch = traverse_pipeline_latches[0];
    if (!latch.valid) {
        return;
    }
    traverse_idle = false;
    // Compute address.
    addr_t dram_address = latch.dram_address;
    if (dram->WillAcceptTransaction(dram_address, false)) {
        dram->AddTransaction(dram_address, false);
        // Create a syndrome for the detector array.
        uint syndrome_size = decoder->circuit.count_detectors()
                                + decoder->circuit.count_observables();
        std::vector<uint8_t> syndrome(syndrome_size, 0);
        for (uint d : latch.detector_array) {
            syndrome[d] = 1;
        }
        dram_await_array.push_back(dram_address);
        addr_to_syndrome[dram_address] = syndrome;
#ifdef GSIM_DEBUG
        std::cout << "\t[Gulliver] Requesting result at address " << 
            std::hex << dram_address << std::dec << " for detectors:";
        for (uint d : latch.detector_array) {
            std::cout << " " << d;
        }
        std::cout << "\n";
#endif
        latch.is_stalled = false;
    } else {
        latch.is_stalled = true;
    }
}

void
GulliverSimulator::tick_traverse_read() {
    TraversePipelineLatch& next = traverse_pipeline_latches[0];
    // Check if we are done.
    if (current_qubit >= decoding_graph.count_detectors()-1) {
        next.valid = false;
        return;
    }
    traverse_idle = false;
    if (next.is_stalled) {
        return;
    }
    DecodingGraph::Vertex currv = decoding_graph.get_vertex(current_qubit);

    std::map<uint, uint> search_region;
    uint search_index = 0;

    std::map<DecodingGraph::Vertex, uint> distance;
    distance[currv] = 0;
    // Perform BFS to determine search region.
    // TODO properly simulate this step.
    std::set<uint> visited;
    std::deque<DecodingGraph::Vertex> bfs_queue{currv};
    while (!bfs_queue.empty()) {
        DecodingGraph::Vertex v = bfs_queue.front();
        bfs_queue.pop_front();
        if (visited.count(v.detector) 
            || distance[v] > search_width) 
        {
            continue;
        }
        search_region[search_index++] = v.detector;
        for (auto w : decoding_graph.adjacency_list(v)) {
            if (w.detector != BOUNDARY_INDEX
                && w.detector >= currv.detector) 
            {
                bfs_queue.push_back(w);
                if (!distance.count(w) || distance[v] + 1 < distance[w]) {
                    distance[w] = distance[v]+1;
                }
            }
        }
        visited.insert(v.detector);
    }
#ifdef GSIM_DEBUG
//  std::cout << "\tSearch region has " << search_region.size() 
//              << " detectors.\n";
#endif
    // Check if the region is all zeros.
    std::vector<uint> detector_array;
    addr_t dram_address = 0x0;
    uint num_nonzeros = 0;
    for (auto pair : search_region) {
        uint index = pair.first;
        uint detector = pair.second;
        if (syndrome_register[detector]) {
            detector_array.push_back(detector); 
#ifdef GSIM_DEBUG
            if (index > (1 << 6)) {
                std::cout << "\tindex of detector in search region "
                    << "is larger than expected: " 
                    << index << "\n";
            }
#endif
            dram_address |= (index & 0x3f) << (num_nonzeros * 6);
            num_nonzeros++;
        }
    }
    dram_address |= current_qubit << 24;

    if (num_nonzeros > max_active_syndromes_in_box) {
        max_active_syndromes_in_box = num_nonzeros;
    }

    if (num_nonzeros > 0 
        && ((num_nonzeros & 0x1) == (hamming_weight_register & 0x1)
            || hamming_weight_register & 0x1)) 
    {
#ifdef GSIM_DEBUG
        if (num_nonzeros > 8) {
            std::cout << "\tSearch region has more than 8 nonzeros: "
                << num_nonzeros << " nonzeros\n";
        }
#endif
        // Update latch.
        next.valid = true;
        next.dram_address = dram_address;
        next.detector_array = detector_array; 
    } else {
        next.valid = false;
    }
    
    current_qubit++;
}

void
GulliverSimulator::update_state() {
    switch (state) {
    case State::traverse:
        state = State::combine;
#ifdef GSIM_DEBUG
        std::cout << "COMBINE\n";
#endif
        break;
    case State::combine:
        state = State::idle;
#ifdef GSIM_DEBUG
        std::cout << "IDLE\n";
#endif
        break;
    default:
        break;
    }
}

void
GulliverSimulator::clear() {
    syndrome_register.clear();
    hamming_weight_register = 0;
    dram_await_array.clear();

    state = State::traverse;
#ifdef GSIM_DEBUG
    std::cout << "TRAVERSE\n";
#endif

    current_qubit = 0;
    for (TraversePipelineLatch& latch : traverse_pipeline_latches) {
        latch.valid = 0;
    }
    cumulative_matching_register.clear();
    matchings.clear();

    addr_to_syndrome.clear();
}

}   // gulliver
}   // qrc
