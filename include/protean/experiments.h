/*
 *  author: Suhas Vittal
 *  date:   2 March 2024
 * */

#ifndef PROTEAN_EXPERIMENTS_h
#define PROTEAN_EXPERIMENTS_h

#include <qontra/graph.h>
#include <qontra/tables.h>

#include <string>

#include "protean/defs.h"

namespace protean {

typedef qgr::Graph<qgr::base::vertex_t, qgr::base::edge_t> CouplingGraph;

uptr<CouplingGraph> read_coupling_graph(std::string);
void make_error_and_timing_from_coupling_graph(std::string, qon::ErrorTable&, qon::TimeTable&);

}   // protean

#endif  // PROTEAN_EXPERIMENTS_h
