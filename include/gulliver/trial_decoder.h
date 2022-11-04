/*
 *  author: Suhas Vittal
 *  date:   1 November 2022
 * */

#ifndef GULLIVER_TRIAL_DECODER_h
#define GULLIVER_TRIAL_DECODER_h

#include "defs.h"
#include "decoding_graph.h"
#include "decoder.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace qrc {
namespace gulliver {

#define GTR_DEBUG

class TrialDecoder : public Decoder {
public:
    TrialDecoder(const stim::Circuit&, uint detectors_per_round,
            const std::filesystem::path& table_directory);
    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;
protected:
    void load_from_file(std::ifstream&);

    typedef std::pair<std::set<DecodingGraph::Edge>, fp_t> ErrorEvent;

    std::map<std::vector<uint8_t>, std::vector<ErrorEvent>> init_error_table;
    std::map<std::vector<uint8_t>, std::vector<ErrorEvent>> middle_error_table;
    std::map<std::vector<uint8_t>, std::vector<ErrorEvent>> final_error_table;
    uint detectors_per_round;
private:
    uint max_candidates;
};

}  // gulliver
}  // qrc

#endif
