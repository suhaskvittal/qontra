/*
 *  author: Suhas Vittal
 *  date:   2 March 2024
 * */

#ifndef QONTRA_DEPENDENCY_GRAPH_h
#define QONTRA_DEPENDENCY_GRAPH_h

#include "qontra/ext/qes.h"
#include "qontra/graph.h"

namespace qontra {
namespace graph {

namespace inst {

struct vertex_t : base::vertex_t {
    qes::Instruction<> instruction;
    size_t meas_order;
};

}   // inst

class DependencyGraph : public Graph<inst::vertex_t, base::edge_t> {
public:
    DependencyGraph(const qes::Program<>&);
    DependencyGraph(DependencyGraph&&) = default;

    void push_back_instruction(const qes::Instruction<>&);

    void contract(sptr<vertex_t>, sptr<vertex_t>);

    qes::Program<> to_program(void) const;
    void remove_redundant_gates(void);

    std::vector<qes::Program<>> get_layers(void) const;

    size_t                  get_depth(void) const;
    sptr<inst::vertex_t>    get_root(void) const;
protected:
    bool update_state(void) override;
private:
    sptr<inst::vertex_t> get_ancestor_of(int64_t);

    std::map<uint64_t, sptr<inst::vertex_t>> recent_ancestor_map;
    sptr<inst::vertex_t> root;
    size_t inst_ctr;
    size_t meas_inst_ctr;

    size_t depth;
};

}   // graph
}   // qontra

#endif  // QONTRA_DEPENDENCY_GRAPH_h
