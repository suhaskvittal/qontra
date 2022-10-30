/*
 *  author: Suhas Vittal
 *  date:   27 October 2022
 * */

#include "gulliver.h"

namespace qrc {

Gulliver::Gulliver(const stim::Circuit circuit, const GulliverParams& params)
    :MWPMDecoder(circuit),
    n_total_accesses(0),
    max_active_syndromes_in_box(0),
    max_hamming_weight(0),
    simulator(nullptr),
    main_clock_frequency(params.main_clock_frequency),
    dram_clock_frequency(params.dram_clock_frequency),
    base_decoder(nullptr),
    dram(nullptr),
    memory_event_table(nullptr)
{
    memory_event_table = new std::map<addr_t, bool>();
    auto cb = [this](addr_t x) 
    {
        this->memory_event_table->insert_or_assign(x, true);
        this->n_total_accesses++;
    };

    dram = new dramsim3::MemorySystem(params.dram_config_file,
                                        params.log_output_directory,
                                        cb, cb);
    base_decoder = new MWPMDecoder(circuit);
    simulator = new gulliver::GulliverSimulator(
                    dram,
                    memory_event_table,
                    graph,
                    base_decoder,
                    params.search_width
                );
}

Gulliver::~Gulliver() {
    dram->PrintStats();

    delete simulator;
    delete base_decoder;
    delete dram;
    delete memory_event_table;
}

std::string
Gulliver::name() {
    return "Gulliver";
}

bool
Gulliver::is_software() {
    return false;
}

DecoderShotResult
Gulliver::decode_error(const std::vector<uint8_t>& syndrome) {
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();

    uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
    if (hw > max_hamming_weight) {
        max_hamming_weight = hw;
    }
    if (hw == 0) {
        return (DecoderShotResult) {
            0.0,
            0.0,
            false,
            std::vector<uint8_t>(),
            std::map<uint, uint>()
        };
    } 
    std::vector<uint8_t> syndrome_dpart(n_detectors, 0);
    for (uint32_t i = 0; i < n_detectors; i++) {
        syndrome_dpart[i] = syndrome[i];
    }
    // Run simulator.
#ifdef GSIM_DEBUG
    std::cout << "==============================================\n";
#endif
    simulator->reset_stats();
    simulator->load_syndrome(syndrome_dpart);

    uint64_t n_cycles = 0;

    fp_t min_clock_frequency = main_clock_frequency < dram_clock_frequency
                                ? main_clock_frequency : dram_clock_frequency;
    uint32_t main_tpc = (uint32_t) (main_clock_frequency / min_clock_frequency);
    uint32_t dram_tpc = (uint32_t) (dram_clock_frequency / min_clock_frequency);

    while (!simulator->is_idle()) {
        for (uint32_t i = 0; i < dram_tpc; i++) {
            dram->ClockTick();    
        }
        for (uint32_t i = 0; i < main_tpc; i++) {
            simulator->tick();
        } 
        n_cycles++;
    }

    auto matching = simulator->get_matching();
    std::vector<uint8_t> correction = get_correction_from_matching(matching);

    fp_t time_taken = n_cycles / min_clock_frequency * 1e9;
    if (simulator->max_active_syndromes_in_box > max_active_syndromes_in_box) {
        max_active_syndromes_in_box = simulator->max_active_syndromes_in_box;
    }

    bool is_error = 
        is_logical_error(correction, syndrome, n_detectors, n_observables);
#ifdef GSIM_DEBUG
    if (is_error) {
        DecoderShotResult mwpm_res = MWPMDecoder::decode_error(syndrome);
        if (!mwpm_res.is_logical_error) {
            std::cout << "Detected correctable error. Dumping matchings.\n";
            std::cout << "[Gulliver]\n";
            for (auto pair : matching) {
                std::cout << "\t" << pair.first << " --> " << pair.second << "\n";
            }
            std::cout << "[MWPM]\n";
            for (auto pair : mwpm_res.matching) {
                std::cout << "\t" << pair.first << " --> " << pair.second
                    << "\terror chain size: " 
                    << graph.get_chain_length(pair.first, pair.second) 
                    << "(" << path_table[pair].distance
                    << ")\n";
            }
        }
    }
#endif

    return (DecoderShotResult) {
        time_taken,
        0.0,
        is_error,
        correction,
        matching
    };
}

} // qrc
