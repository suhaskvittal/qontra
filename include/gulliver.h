/*
 *  author: Suhas Vittal
 *  date:   25 October 2022
 * */

#ifndef GULLIVER_h
#define GULLIVER_h

#include "defs.h"
#include "mwpm_decoder.h"
#include "gulliver/simulator.h"

#include <algorithm>
#include <map>
#include <string>

namespace qrc {

struct GulliverParams {
    uint search_width;  // Should be (d+1)/2

    std::string dram_config_file;
    std::string log_output_directory;

    fp_t main_clock_frequency;
    fp_t dram_clock_frequency;
};

class Gulliver : public MWPMDecoder {
public:
    Gulliver(const stim::Circuit, const GulliverParams&);
    ~Gulliver();

    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;

    // Statistics
    uint32_t n_total_accesses;
    uint32_t max_active_syndromes_in_box;
    uint max_hamming_weight;

    gulliver::GulliverSimulator * simulator;
private:
    fp_t main_clock_frequency;
    fp_t dram_clock_frequency;

    MWPMDecoder * base_decoder;
    dramsim3::MemorySystem * dram;
    std::map<addr_t, bool> * memory_event_table;
};

} // qrc

#endif
