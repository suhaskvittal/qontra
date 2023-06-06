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

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace qontra {
namespace protean {

namespace compiler {
    typedef std::function<bool(Processor3D&)> constraint_t;   // Constraint type: outputs
                                                              // true if the graph obeys
                                                              // the specified constraint.
    typedef std::function<fp_t(Processor3D&)> cost_t;         // Returns a float scoring
                                                              // the graph. The compiler
                                                              // will try and minimize this
                                                              // score.
}   // compiler

class Compiler {
public:
    typedef struct {
    } params_t;

    typedef struct {
        TannerGraph                     curr_spec;
        Processor3D                     arch;
        std::vector<qc::Instruction>    schedule;
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

    Compiler(const std::vector<compiler::constraint_t>& con, compiler::cost_t obj)
        :constraints(con), 
        objective(obj),
        max_induced_check_weight(std::numeric_limits<uint>::max()),
        verbosity(0),
        compile_round(0)
    {}

    ir_t    run(const TannerGraph&, bool verbose=true);
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
    //      (8) Linearize -- linearize data connections in the exist architecture. Jump to (2).
    //
    //  We repeat the following until we exit at (4). We assume that the optimization space is
    //  "smooth", so modifications do not chaotically affect the score.
    void    place(ir_t&);
    void    reduce(ir_t&);
    void    merge(ir_t&);
    void    schedule(ir_t&);
    void    score(ir_t&);
    bool    induce(ir_t&);
    bool    sparsen(ir_t&);
    void    linearize(ir_t&);

    const std::vector<compiler::constraint_t>   constraints;
    const compiler::cost_t                      objective;

    uint max_induced_check_weight;  // We will not create induced predecessors for checks above
                                    // this weight.
    uint verbosity;

    uint compile_round;
    bool called_sparsen;
};

void    print_connectivity(Processor3D&);
void    print_schedule(const std::vector<qc::Instruction>&);

}
}   // qontra

#endif  // PROTEAN_COMPILER_h
