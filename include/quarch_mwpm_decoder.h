/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef QUARCH_MWPM_DECODER_h
#define QUARCH_MWPM_DECODER_h

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

#include <stim.h>

#include <PerfectMatching.h>

#include "quarch_decoder.h"
#include "quarch_lilliput.h"

#include <vector>
#include <map>
#include <chrono>
#include <set>
#include <utility>

#include <time.h>

#define MWPM_DECODER_NAME "MWPMDecoder"
#define MAX_PRACTICAL_DISTANCE  1000.0

#if QFP_SIZE==0
#define MWPM_INTEGER_SCALE 10000.0
#elif QFP_SIZE==1
#define MWPM_INTEGER_SCALE 1000.0
#elif QFP_SIZE==2
#define MWPM_INTEGER_SCALE 10.0
#endif

struct DijkstraResult {
    std::vector<uint> path;
    fp_t distance;
};

typedef qfp_t  wgt_t;

class MWPMDecoder : public Decoder {
public:
    MWPMDecoder(const stim::Circuit&);

    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;

    uint64_t sram_cost(void) override;
protected:
    std::map<std::pair<uint, uint>, DijkstraResult> path_table;
private:
    void update_path_table(uint src, uint dst, 
            const std::vector<fp_t>& distances,
            const std::vector<uint>& predecessors);
};

#endif
