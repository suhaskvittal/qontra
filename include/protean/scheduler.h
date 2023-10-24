/*
 *  author: Suhas Vittal
 *  date:   23 October 2023
 * */

#ifndef PROTEAN_SCHEDULER_h
#define PROTEAN_SCHEDULER_h

#include "defs.h"
#include "graph/tanner_graph.h"

#include <ilcplex/ilocplex.h>

#include <deque>
#include <map>
#include <set>
#include <vector>

namespace qontra {
namespace protean {

#define __QUBIT_LIST    std::vector<graph::tanner::vertex_t*>
#define __CHECK_DEFS    std::map<graph::tanner::vertex_t*, __QUBIT_LIST>
#define __FLAG_DEFS     std::map<graph::tanner::vertex_t*, graph::tanner::vertex_t*>

struct css_code_data_t {
    std::vector<uint>   data_qubits;
    std::vector<uint>   parity_qubits;
    std::vector<uint>   flag_qubits;
    
    // Each entry of check schedules should be an array of integers such that the i-th
    // entry is the qubit that will have a CNOT with the parity qubit at time step i.
    // 
    // This should be -1 if no such qubit exists.
    std::map<uint, std::vector<int32_t>>    check_schedules;
    
    std::vector<uint>   xparity_list;
    std::vector<uint>   zparity_list;

    std::vector<std::vector<uint>>  x_obs_list;
    std::vector<std::vector<uint>>  z_obs_list;

    // flag_usage[u][v] corresponds to the flag qubit used by data qubit v during check u.
    TwoLevelMap<uint, uint, uint>   flag_usage;
};

// Each stabilizer is defined as an array of pauli_ops.
// We choose such a definition to allow this code to be
// extensible to generic stabilizer codes.
enum class pauli { x, y, z };
typedef std::pair<pauli, uint>  pauli_op_t;
typedef std::vector<pauli_op_t> stabilizer_t;
typedef std::vector<uint>       support_t;

// We define a general stabilizer graph which connects
// two stabilizers that have anticommuting pauli operators
// in their support (note that by definition, the stabilizers
// themselves commute -- but an operator on a single qubit
// may not).
struct stab_vertex_t : graph::base::vertex_t {
    tanner::vertex_t* check;
    stabilizer_t stabilizer;
    support_t support_t;
    std::map<uint, pauli> qubit_to_pauli;

    // Schedule maps qubit to time. This is easier for processing later.
    std::map<uint, uint> sch_qubit_to_time;
    std::map<uint, uint> sch_time_to_qubit;
};

struct stab_edge_t : graph::base::edge_t {
};

typedef graph::Graph<stab_vertex_t, stab_edge_t>    StabilizerGraph;

// Finally, the scheduling function.
css_code_data_t
compute_schedule_from_tanner_graph(graph::TannerGraph&);

// Helper functions:
IloModel
construct_scheduling_program(stab_vertex_t*,
                                StabilizerGraph&, 
                                int max_stabilizer_weight,
                                std::map<uint, IloIntVar>& qubit_to_variable);

void    finalize_cplex(void);

}   // protean
}   // qontra

#endif  // PROTEAN_SCHEDULER_h
