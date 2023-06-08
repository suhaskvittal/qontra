/*
 *  author: Suhas Vittal
 *  date:   31 May 2023
 * */

#ifndef PROTEAN_COMPILER_h
#define PROTEAN_COMPILER_h

#include "defs.h"
#include "instruction.h"
#include "protean/proc3d.h"
#include "protean/tanner_graph.h"

#include <lemon/list_graph.h>
#include <lemon/matching.h>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace qontra {
namespace protean {

namespace compiler {
    typedef struct {
        TannerGraph                     curr_spec;
        Processor3D                     arch;
        schedule_t<qc::Instruction>     schedule;
        fp_t                            score;
        bool                            valid;
        // Data structures:
        std::map<tanner::vertex_t*, proc3d::vertex_t*>              role_to_qubit;
        std::map<proc3d::vertex_t*, std::vector<tanner::vertex_t*>> qubit_to_roles;
                                                        // During the "reduce"  pass, a qubit
                                                        // may take on multiple roles due to a
                                                        // contraction. For example, a parity
                                                        // check qubit may need to take on the
                                                        // role of a gauge qubit later.
        std::set<proc3d::vertex_t*> is_gauge_only;  // Keep track of pure gauge qubits for
                                                    // operations like reduce.
    } ir_t;

    typedef std::function<bool(ir_t*)> constraint_t;   // Constraint type: outputs
                                                              // true if the graph obeys
                                                              // the specified constraint.
    typedef std::function<fp_t(ir_t*)> cost_t;         // Returns a float scoring
                                                              // the graph. The compiler
                                                              // will try and minimize this
                                                              // score.
}   // compiler

class Compiler {
public:
    typedef struct {
    } params_t;

    Compiler(const std::vector<compiler::constraint_t>& con, compiler::cost_t obj)
        :constraints(con), 
        objective(obj),
        max_induced_check_weight(std::numeric_limits<uint>::max()),
        verbosity(0),
        compile_round(0)
    {}

    compiler::ir_t* run(const TannerGraph&, bool verbose=true);
private:
    // Compiler passes:
    //  (1) Place       -- creates a architectural description for the current Tanner graph.
    //  (2) Merge       -- merges parity checks that act on the same qubits (i.e. color code)
    //  (3) Reduce      -- removes redundant qubits from the architecture through contraction.
    //  (4) Schedule    -- schedules the operations to compute each check.
    //  (5) Score       -- scores the architecture: if the score worsens, exit.
    //  (6) Induce      -- induces predecessors onto the Tanner graph.
    //  If Induce fails:
    //      (7) Sparsen -- if a gauge/parity qubit has too many connections, then split them
    //                      into groups of 2 by using intermediate gauge qubits. 
    //  If Sparsen fails:
    //      (8) Linearize -- linearize data connections in the exist architecture. Jump to (3).
    //
    //  We repeat the following until we exit at (4). We assume that the optimization space is
    //  "smooth", so modifications do not chaotically affect the score.
    void    place(compiler::ir_t*);
    void    reduce(compiler::ir_t*);
    void    merge(compiler::ir_t*);
    void    schedule(compiler::ir_t*);
    void    score(compiler::ir_t*);
    bool    induce(compiler::ir_t*);
    bool    sparsen(compiler::ir_t*);
    void    linearize(compiler::ir_t*);

    const std::vector<compiler::constraint_t>   constraints;
    const compiler::cost_t                      objective;

    uint max_induced_check_weight;  // We will not create induced predecessors for checks above
                                    // this weight.
    uint verbosity;

    uint compile_round;
    bool called_sparsen;
};

void    print_connectivity(Processor3D&);
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
