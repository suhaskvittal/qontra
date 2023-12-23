/*
 *  author: Suhas Vittal
 *  date:   27 November 2023
 * */

#ifndef QONTRA_STIM_EXTENSIONS
#define QONTRA_STIM_EXTENSIONS

#include <stim.h>

#include <map>
#include <set>

namespace qontra {

// We provide an extension of Stim's Circuit class. The DetailedStimCircuit
// class can be used in any function accepting stim::Circuits, but has
// extra metadata that can be used by decoders and such.
//
// Stim's circuits cannot hold much information beyond error rates
// and operations, so we contain extra information in here.
struct DetailedStimCircuit : public stim::Circuit {
    DetailedStimCircuit() :stim::Circuit() {}
    DetailedStimCircuit(const stim::Circuit& other)
        :stim::Circuit(other)
    {}
    DetailedStimCircuit(const stim::Circuit&& other)
        :stim::Circuit(std::move(other))
    {}
    DetailedStimCircuit(const DetailedStimCircuit& other)
        :stim::Circuit(other),
        detection_event_to_color(other.detection_event_to_color),
        flag_detection_events(other.flag_detection_events),
        flag_edge_table(other.flag_edge_table)
    {}

    DetailedStimCircuit(const DetailedStimCircuit&& other)
        :stim::Circuit(other),
        detection_event_to_color(std::move(other.detection_event_to_color)),
        flag_detection_events(std::move(other.flag_detection_events)),
        flag_edge_table(other.flag_edge_table)
    {}

    typedef std::tuple<uint, uint, uint>    flag_edge_t;

    std::map<uint, int> detection_event_to_color;
    std::set<uint>      flag_detection_events;

    std::map<uint, flag_edge_t> flag_edge_table;
};

}   // qontra

#endif  // QONTRA_STIM_EXTENSIONS
