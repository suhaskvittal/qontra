/*
 *  author: Suhas Vittal
 *  date:   23 October 2023
 * */

#ifndef PROTEAN_SCHEDULER_h
#define PROTEAN_SCHEDULER_h

#include "defs.h"
#include "experiments.h"
#include "graph/tanner_graph.h"
#include "instruction.h"
#include "linprog.h"
#include "protean/representation.h"
#include "protean/utils.h"
#include "sim/enumerator.h"

#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <assert.h>
#include <math.h>
#include <string.h>

namespace qontra {
namespace protean {

// Each stabilizer is defined as an array of pauli_ops.
// We choose such a definition to allow this code to be
// extensible to generic stabilizer codes.

#define pauli   enumerator::error

typedef std::pair<pauli, uint>  pauli_op_t;
typedef std::vector<pauli_op_t> mq_pauli_op_t;
typedef mq_pauli_op_t           stabilizer_t;
typedef std::vector<uint>       support_t;

mq_pauli_op_t
mul(const mq_pauli_op_t&, const mq_pauli_op_t&);

int
count_anticommutations(const mq_pauli_op_t&, const mq_pauli_op_t&);

std::vector<stabilizer_t>   
get_stabilizers(css_code_data_t); 

// We define a general stabilizer graph which connects
// two stabilizers that have anticommuting pauli operators
// in their support (note that by definition, the stabilizers
// themselves commute -- but an operator on a single qubit
// may not).

struct stab_vertex_t : graph::base::vertex_t {
    graph::tanner::vertex_t* check;
    stabilizer_t stabilizer;
    support_t support;
    std::map<uint, pauli> qubit_to_pauli;

    // Schedule maps qubit to time. This is easier for processing later.
    std::map<uint, uint> sch_qubit_to_time;
    std::map<uint, uint> sch_time_to_qubit;
};

struct stab_edge_t : graph::base::edge_t {
};

typedef graph::Graph<stab_vertex_t, stab_edge_t>    StabilizerGraph;

// compute_schedule_from_tanner_graph computes a functionally valid schedule
// with depth that is theoretically guaranteed to be at worst 1.5x longer than
// the shortest depth schedule.
//
// This schedule is computed by traversing across the stabilizers and scheduling
// each stabilizer with respect to prior scheduled stabilizers.
//
// The "start" parameter is optional, but we found that different choices of
// starting stabilizer can yield different results. As the function is fast,
// it is indeed possible to iterate over all possible starting points and choose
// the best result.
css_code_data_t
compute_schedule_from_tanner_graph(graph::TannerGraph&, int start=0);

css_code_data_t
make_fault_tolerant(css_code_data_t);

// Helper functions:
LPManager<uint>*
construct_scheduling_program(stab_vertex_t*,
                                StabilizerGraph&, 
                                int max_stabilizer_weight);

}   // protean
}   // qontra

#endif  // PROTEAN_SCHEDULER_h
