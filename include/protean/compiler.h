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
    typedef std::function<bool(const Processor3D&)> constraint_t;   // Constraint type: outputs
                                                                    // true if the graph obeys
                                                                    // the specified constraint.
    typedef std::function<fp_t(const Processor3D&)> cost_t;         // Returns a float scoring
                                                                    // the graph. The compiler
                                                                    // will try and minimize this
                                                                    // score.
}   // compiler

class Compiler {
public:
    typedef struct {
    } params_t;

    typedef struct { 
        fp_t                            score; 
        bool                            valid; 
        Processor3D                     arch;
        std::vector<qc::Instruction>    schedule;
    } result_t;

    Compiler(const std::vector<compiler::constraint_t>& con, compiler::cost_t obj)
        :constraints(con), 
        objective(obj),
        max_induced_check_weight(std::numeric_limits<uint>::max())
    {}

    result_t    run(const TannerGraph&);
private:
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
    } ir_t;

    // Compiler passes:
    //  (1) Place       -- creates a architectural description for the current Tanner graph.
    //  (2) Reduce      -- removes redundant qubits from the architecture through contraction.
    //  (3) Schedule    -- schedules the operations to compute each check.
    //  (4) Score       -- scores the architecture: if the score worsens, exit.
    //  (5) Induce      -- induces predecessors onto the Tanner graph.
    //  If Induce fails:
    //      (6) Sparsen -- if a gauge/parity qubit has too many connections, then split them
    //                      into groups of 2 by using intermediate gauge qubits. 
    //  If Sparsen fails:
    //      (7) Linearize -- linearize data connections in the exist architecture. Jump to (2).
    //
    //  We repeat the following until we exit at (4). We assume that the optimization space is
    //  "smooth", so modifications do not chaotically affect the score.
    void    place(ir_t&);
    void    reduce(ir_t&);
    void    schedule(ir_t&);
    void    score(ir_t&);
    bool    induce(ir_t&);
    bool    sparsen(ir_t&);
    void    linearize(ir_t&);

    const std::vector<compiler::constraint_t>   constraints;
    const compiler::cost_t                      objective;

    uint max_induced_check_weight;  // We will not create induced predecessors for checks above
                                    // this weight.
};

}
}   // qontra

#endif  // PROTEAN_COMPILER_h
