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

#include <functional>
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
    Compiler(const std::vector<compiler::constraint_t>& con, compiler::cost_t obj)
        :constraints(con), objective(obj)
    {}

    typedef struct { 
        fp_t                            score; 
        bool                            valid; 
        Processor3D                     arch;
        std::vector<qc::Instruction>    schedule;
    } result_t;

    result_t    run(const TannerGraph&);
private:
    typedef struct {
        TannerGraph                     curr_spec;
        Processor3D                     arch;
        std::vector<qc::Instruction>    schedule;

        fp_t                            score;
        bool                            valid;
    } ir_t;

    // Compiler passes:
    //  (1) Place       -- creates a architectural description for the current Tanner graph.
    //  (2) Loop:
    //      (a) Reduce      -- removes redundant qubits from the architecture through contraction.
    //      (b) Schedule    -- schedules the operations to compute each check.
    //  (4) Score       -- scores the architecture: if the score worsens, exit.
    //  (5) Induce      -- induces predecessors onto the Tanner graph.
    //
    //  We repeat the following until we exit at (4). We assume that the optimization space is
    //  "smooth", so modifications do not chaotically affect the score.
    void    place(ir_t&);
    void    reduce(ir_t&);
    void    schedule(ir_t&);
    void    score(ir_t&);
    void    induce(ir_t&);

    const std::vector<compiler::constraint_t>   constraints;
    const compiler::cost_t                      objective;

    result_t    best_result;
};
    

}
}   // qontra

#endif  // PROTEAN_COMPILER_h
