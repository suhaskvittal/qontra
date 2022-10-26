/*
 *  author: Suhas Vittal
 *  date:   22 September 2022
 * */

#include "hyperion/cache.h"

namespace qrc {
namespace hyperion {

static QubitCache::SuperTag DEFAULT_SUPERTAG = {0, 0, false};
static QubitCache::Tag DEFAULT_TAG = {0, 0, 0, false, true, 0}; 

QubitCache::QubitCache(uint n_supertags, uint n_sets, uint n_lines, uint min_detector)
    :load_queue(),
    replacement_queue(),
    supertag_store(n_supertags, DEFAULT_SUPERTAG),
    sets(n_sets, std::vector<Tag>(n_lines/n_sets, DEFAULT_TAG)),
    result_table(),
    min_detector(min_detector)
{}

void
QubitCache::tick() {
    // Perform a replacement.
    if (!replacement_queue.empty()) {
        // Note that we don't need the data.
        // We put the tag into the cache, and we load
        // the data only during the ccomp state of the
        // processor.
        //
        // Any hits on the cache will be an effective 
        // miss since the data is marked as incomplete.
        auto new_entry = replacement_queue.front();
        uint new_qubit = std::get<0>(new_entry);
        uint new_detector = std::get<1>(new_entry);
        bool super_replace = std::get<2>(new_entry);

        n_evictions++;
        if (super_replace) {
            n_super_evictions++;
        }
        replace(new_qubit, new_detector, super_replace);
        replacement_queue.pop_front();
    }
    // Add an entry into the cache.
    if (!load_queue.empty()) {
        n_accesses++;

        auto request = load_queue.front();
       
        AccessOutcome outcome = is_hit(request);
        if (outcome != AccessOutcome::hit) {
            n_misses++;
            // Then we have to schedule a replacement.
            bool is_super_miss = (outcome == AccessOutcome::supermiss);
            replacement_queue.push_back(
                    std::make_tuple(request.qubit, request.di, is_super_miss));
            if (request.dj >= min_detector) {
                replacement_queue.push_back(
                        std::make_tuple(request.qubit, request.dj, is_super_miss));
            }
            // Now, the processor must load the data for the replacement to
            // occur.
            result_table[request] = LoadStatus::miss;
        } else {
            result_table[request] = LoadStatus::hit;
        }
    }
}

QubitCache::LoadStatus
QubitCache::access(uint qubit, uint di, uint dj) {
    if ((di < min_detector || di == BOUNDARY_INDEX) 
            && (dj < min_detector || dj == BOUNDARY_INDEX)) 
    {
        return LoadStatus::miss;
    }

    if (di != BOUNDARY_INDEX && di < dj) {  // Also covers when di < min_detector.
        // Enforce ordering (di > dj).
        return access(qubit, dj, di);
    }

    Request k = {qubit, di, dj};
    if (result_table.count(k)) {  // This is an old access.
        LoadStatus s = result_table[k];
        result_table.erase(k);
        return s;
    } else {
        // This requires an access. 
        load_queue.push_back(k);
        result_table[k] = LoadStatus::unknown;
        return LoadStatus::unknown;
    }
}

void
QubitCache::complete(uint qubit, uint detector) {
    uint si = qubit % sets.size();
    for (Tag& t : sets[si]) {
        if (t.qubit == qubit && t.detector == detector && t.valid) {
            t.incomplete = false;
            break;
        }
    }
}

void
QubitCache::replace(uint qubit, uint detector, bool super_replace) {
    // If this is a super replace, find a victim from the supertag store.
    // We want to choose a victim whose logical qubit maps to the same
    // set as our new logical qubit (otherwise, what's the point?)
    uint si = qubit % sets.size();
    if (super_replace) {
        uint victim = 0;
        uint32_t min_frequency = (uint32_t) -1;

        for (uint i = 0; i < supertag_store.size(); i++) {
            const SuperTag& st = supertag_store[i];
            if (!st.valid) {
                victim = i;
                break;
            } else if (si == (st.qubit % sets.size())
                    && st.frequency < min_frequency) 
            {
                victim = i;
                min_frequency = st.frequency;
            }
        }
        
        // First, invalid all tags corresponding to the victim supertag.
        SuperTag& victim_tag = supertag_store[victim];
        uint i = victim_tag.qubit % sets.size();
        for (Tag& t : sets[i]) {
            if (t.qubit == victim_tag.qubit) {
                t.valid = false;
                t.incomplete = true;
            }
        }
    }
    // Replace a tag from the set.
    uint victim = 0;
    uint32_t min_frequency = (uint32_t) -1;

    for (uint i = 0; i < sets[si].size(); i++) {
        Tag& t = sets[si][i]; 
        if (t.rprot_cycles_left--) {
            continue;
        }
        if (!t.valid) {
            victim = i;
            break;
        } else if (t.frequency < min_frequency) {
            victim = i;
            min_frequency = t.frequency;
        }
    }

    sets[si][victim] = {
        qubit,
        detector,
        true,
        true,
        1,
        RPROT_CYCLES
    };
    completion_queue.push_back(detector);
}

QubitCache::AccessOutcome
QubitCache::is_hit(const Request& r) {
    n_accesses++;
    
    bool is_miss;
    // Check if we can match the super tag.
    is_miss = true;
    for (SuperTag& st : supertag_store) {
        if (st.valid && st.qubit == r.qubit) {
            st.frequency++;
            is_miss = false;
            break;
        }
    }

    if (is_miss) {
        return AccessOutcome::supermiss;
    }

    is_miss = true;
    // Check if we can match the tag.
    uint si = r.qubit % sets.size();
    for (Tag& t : sets[si]) {
        if (t.qubit == r.qubit 
                && (t.detector == r.di || t.detector == r.dj)
                && t.valid)
        {
            t.frequency++;
            if (!t.incomplete) {
                is_miss = false;
            }
        }
    }

    if (is_miss) {
        return AccessOutcome::miss;
    }
    // Otherwise, this was a hit.
    return AccessOutcome::hit;
}

}
}
