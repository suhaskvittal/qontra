/* author: Suhas Vittal
 *  date:   31 May 2023
 * */

#ifndef PROTEAN_COMPILER_h
#define PROTEAN_COMPILER_h

#include "defs.h"
#include "graph/algorithms/distance.h"
#include "graph/algorithms/mis.h"
#include "graph/algorithms/search.h"
#include "graph/dependence_graph.h"
#include "graph/graph.h"
#include "graph/tanner_graph.h"
#include "instruction.h"
#include "protean/proc3d.h"

#include <lemon/list_graph.h>
#include <lemon/matching.h>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <vector>
#include <utility>

namespace qontra {
namespace protean {

namespace compiler {
    struct ir_t {
        ~ir_t(void) {
            delete arch;
            delete dependency_graph;
        }

        graph::TannerGraph*             curr_spec;
        Processor3D*                    arch;
        schedule_t<qc::Instruction>     schedule;
        fp_t                            score;
        bool                            valid;
        // Data structures:
        std::map<graph::tanner::vertex_t*, proc3d::vertex_t*>
                                                role_to_qubit;
        std::map<proc3d::vertex_t*, std::vector<graph::tanner::vertex_t*>>
                                                qubit_to_roles;
                                                        // During the "reduce"  pass, a qubit
                                                        // may take on multiple roles due to a
                                                        // contraction. For example, a parity
                                                        // check qubit may need to take on the
                                                        // role of a gauge qubit later.
        std::set<proc3d::vertex_t*> is_gauge_only;  // Keep track of pure gauge qubits for
                                                    // operations like reduce.

        std::map<tanner::vertex_t*, schedule_t<qc::Instruction>>    check_to_impl;  
                                                            // Contains micro-schedules which
                                                            // implement each parity check.
        std::set<std::pair<tanner::vertex_t*, tanner::vertex_t*>>   conflicting_checks;
                                                            // A set of checks that cannot be
                                                            // scheduled concurrently.
        graph::DependenceGraph<qc::Instruction>*                    dependency_graph;
                                                            // Defines the dependence relation
                                                            // for the instrucitons in the 
                                                            // macro-schedule.
                                            
        std::set<tanner::vertex_t*> sparsen_visited_set;

        bool is_data(proc3d::vertex_t* v) {
            return !qubit_to_roles.count(v);
        }
    };

    typedef std::function<fp_t(ir_t*)> cost_t;  // Returns a float scoring the IR. The compiler
                                                // will try and minimize this score.
    typedef struct {
        uint max_qubits         = std::numeric_limits<uint>::max();
        uint max_connectivity   = std::numeric_limits<uint>::max();
        uint max_thickness      = std::numeric_limits<uint>::max();

        fp_t max_mean_connectivity =    std::numeric_limits<fp_t>::max();
    } constraint_t;
}   // compiler

class Compiler {
public:
    typedef struct {
        bool verbose =  false;

    } params_t;

    Compiler(compiler::constraint_t con, compiler::cost_t obj)
        :constraints(con), 
        objective(obj),
        max_induced_check_weight(std::numeric_limits<uint>::max()),
        compile_round(0),
        params(),
        rng(0)
    {}

    void            set_seed(uint64_t s) { rng.seed(s); }

    compiler::ir_t* run(graph::TannerGraph*);

    params_t    params;
private:
    // These constraint checks return true if there is a violation.
    bool check_connectivity_violation(compiler::ir_t* ir) {
        return ir->arch->get_mean_connectivity() > constraints.max_mean_connectivity
            || ir->arch->get_max_connectivity() > constraints.max_connectivity;
    }

    bool check_size_violation(compiler::ir_t* ir) {
        return ir->arch->get_vertices().size() > constraints.max_qubits;
    }
    // Compiler passes:
    //  (1) Place       -- creates a architectural description for the current Tanner graph.
    //                      Place is not guaranteed to obey constraints. Instead, future passes
    //                      should try and fit the initial placement according to the specified
    //                      constraints.
    //  (2) Unify       -- merges parity checks that act on the same qubits (i.e. color code)
    //  (3) Reduce      -- removes redundant qubits from the architecture through contraction.
    //
    //  (4)
    //  If constraints.max_qubits is violated:
    //      (a) Merge: Find a pair of adjacent non-data qubits and merge them. Choose this pair
    //                  such that their degrees are minimized.
    //      Jump to (3).
    //  If constraints.max_connectivity is violated:
    //      (b) Split: Find the most connected qubit and split it into two. If the qubit is a 
    //                  a data qubit, then create a gauge qubit.
    //      Jump to (3).
    //
    //  (5) Micro-Schedule  -- schedules the operations for each check.
    //  (6) Macro-Schedule  -- schedules the order of computing each check such that depth is
    //                          minimized.
    //  
    //  (6) Score   -- check if IR is valid. If not, jump to 7.
    //  If Observable is defined:
    //      (a) Find the best CNOT order by repeatedly calling the score function.
    //  Else:
    //      (b) Score normally.
    //
    //  (7) Induce      -- induces predecessors onto the Tanner graph.
    //
    //  If Induce fails:
    //      (8) Sparsen -- if a gauge/parity qubit has too many connections, then split them
    //                      into groups of 2 by using intermediate gauge qubits.
    //  If Sparsen fails:
    //      (9) Linearize -- linearize data connections in the exist architecture. Jump to (3).
    //
    //  We repeat the following until we exit at (4). We assume that the optimization space is
    //  "smooth", so modifications do not chaotically affect the score.
    void    place(compiler::ir_t*);
    void    unify(compiler::ir_t*);
    void    reduce(compiler::ir_t*);
    bool    merge(compiler::ir_t*);
    bool    split(compiler::ir_t*);
    void    micro_schedule(compiler::ir_t*);
    void    macro_schedule(compiler::ir_t*);
    void    schedule(compiler::ir_t*);
    void    score(compiler::ir_t*);
    bool    induce(compiler::ir_t*);
    bool    sparsen(compiler::ir_t*);
    void    linearize(compiler::ir_t*);

    const compiler::constraint_t    constraints;
    const compiler::cost_t          objective;

    uint max_induced_check_weight;  // We will not create induced predecessors for checks above
                                    // this weight.
    uint compile_round;
    bool called_sparsen;

    std::mt19937_64 rng;
};

void    print_connectivity(Processor3D*);
void    print_schedule(const schedule_t<qc::Instruction>&);

// write_ir_to_folder dumps an IR to a folder into the following files:
//  (1) spec.txt    (TannerGraph)
//  (2) arch/
//      (a) flat_map.txt    (Processor3D, connections without verticality, not necessarily planar)
//      (b) 3d_map.txt      (Processor3D, connections with verticality, k-planar)
//  (3) labels.txt  (role_to_qubit)
//  (4) schedule.qasm (schedule)
void    write_ir_to_folder(compiler::ir_t*, std::string);

}   // protean
}   // qontra

#endif  // PROTEAN_COMPILER_h
