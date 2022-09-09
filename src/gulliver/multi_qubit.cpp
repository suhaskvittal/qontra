/*
 *  author: Suhas Vittal
 *  date:   9 September 2022
 * */

#include "gulliver/multi_qubit.h"

GulliverMultiQubitSimulator::GulliverMultiQubitSimulator(
        const std::vector<stim::Circuit>& circuits,
        uint n_decoders,
        const GulliverParams& params)
    :n_logical_failures(0),
    dram(nullptr),
    memory_event_table(nullptr),
    simulators(n_decoders),
    circuits(circuits),
    main_clock_frequency(params.main_clock_frequency)
{
    // Initialize memory first.
    memory_event_table = new std::map<addr_t, bool>(); 
    auto cb = [this](addr_t x) 
    {
        this->memory_event_table->insert_or_assign(x, true);
    };
    dram = new dramsim3::MemorySystem(params.dram_config_file,
                                        params.log_output_directory,
                                        cb, cb);

}
