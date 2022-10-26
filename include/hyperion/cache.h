/*
 *  author: Suhas Vittal
 *  date:   7 September 2022
 * */

#ifndef HYPERION_CACHE_h
#define HYPERION_CACHE_h

#include "defs.h"
#include "decoding_graph.h"

#include <deque>
#include <limits>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

namespace qrc {
namespace hyperion {

#define RPROT_CYCLES    8

class QubitCache {
public:
    QubitCache(uint n_supertags, uint n_sets, uint n_lines, uint min_detector);
    
    enum class LoadStatus {hit, miss, unknown};

    struct SuperTag {
        uint qubit; 
        uint32_t frequency;
        bool valid;
    };

    struct Tag {
        uint qubit;
        uint detector;
        bool valid;
        bool incomplete;
        uint32_t frequency;
        uint32_t rprot_cycles_left;
    };

    void tick(void);
    LoadStatus access(uint qubit, uint di, uint dj); // This will either
                                                     // (1) Add the access
                                                     //    to the laod queue,
                                                     // (2) Return the load 
                                                     //    status of the access
                                                     //    if it has been completed.
    void complete(uint qubit, uint detector);

    std::deque<uint> completion_queue;

    // Statistics
    uint32_t n_accesses;
    uint32_t n_misses;
    uint32_t n_super_evictions;
    uint32_t n_evictions;
private:
    enum class AccessOutcome {hit, miss, supermiss};
    struct Request {
        uint qubit;
        uint di;
        uint dj;

        bool operator==(const Request& r) const {
            return qubit == r.qubit
                && di == r.di
                && dj == r.dj;
        }

        bool operator<(const Request& r) const {
            return qubit < r.qubit
                || (qubit == r.qubit && di < r.di)
                || (qubit == r.qubit && di == r.di && dj < r.dj);
        }
    };

    void replace(uint qubit, uint detector, bool super_replace);
    AccessOutcome is_hit(const Request&);

    // Microarchitectural elements
    std::deque<Request> load_queue; // Data in cache is read-only
                                                  // so we don't need a store 
                                                  // queue.
    std::deque<std::tuple<uint, uint, bool>> replacement_queue; // Data in queue is
                                                                // (1) qubit
                                                                // (2) detector
                                                                // (3) level of 
                                                                //     replacement
    std::vector<SuperTag> supertag_store;
    std::vector<std::vector<Tag>> sets;
    // Data
    std::map<Request, LoadStatus> result_table;
    uint min_detector;
};

}
}

#endif
