/*
 *  author: Suhas Vittal
 *  date:   2 March 2024
 * */

#include "qontra/graph/dependency_graph.h"

#include <vtils/utility.h>

#include <deque>
#include <limits>

namespace qontra {
namespace graph {

using namespace inst;

DependencyGraph::DependencyGraph(const qes::Program<>& program)
    :recent_ancestor_map(),
    root(nullptr),
    inst_ctr(0),
    meas_inst_ctr(0)
{
    root = make_and_add_vertex(std::numeric_limits<uint64_t>::max());
    for (size_t i = 0; i < program.size(); i++) {
        const qes::Instruction<>& inst = program.at(i);
        if (!is_gate(inst)) {
            std::cerr << "[ DependencyGraph ] found non-quantum instruction:" 
                << inst << " on line " << i+1 << std::endl;
            exit(1);
        }
        push_back_instruction(inst);
    }
}

qes::Program<>
DependencyGraph::to_program() const {
    std::vector<qes::Program<>> layers = get_layers();
    qes::Program<> program;
    for (auto& blk : layers) {
        vtils::push_back_range(program, blk);
    }
    return program;
}

void
DependencyGraph::remove_redundant_gates() {
    std::deque<sptr<vertex_t>> bfs{root};
    std::map<sptr<vertex_t>, sptr<vertex_t> prev;

    std::set<sptr<vertex_t>> visited;
    while (bfs.size()) {
        sptr<vertex_t> v = bfs.front();
        bfs.pop_front();
        if (visited.count(v)) continue;
        for (sptr<vertex_t> w : get_outgoing(v)) {
            bfs.push_back(w);
            prev[w] = v;
        }
        // Try and merge with the previous vertex.
        if (v != root) {
            sptr<vertex_t> w = prev.at(v);
            if (v->get_name() == w->get_name()
                    && v->get_operands() == w->get_operands())
            {
                contract(v, w);
            }
        }
        visited.insert(v);
    }
}

std::vector<qes::Program<>>
DependencyGraph::get_layers() const {
    std::deque<sptr<vertex_t>> bfs{root};
    std::map<sptr<vertex_t>, int> layer_map;
    layer_map[root] = -1;

    std::set<sptr<vertex_t>> visited;
}

void
DependencyGraph::push_back_instruction(const qes::Instruction<>& inst) {
    size_t k = is_2q_gate(inst) ? 2 : 1;

    if (inst.size() > k) {
        // Split the instruction into multiple subinstructions.
        for (size_t i = 0; i < inst.get_number_of_operands(); i += k) {
            std::vector<int64_t> operands;
            for (size_t j = i; j < i+k; j++) {
                operands.push_back(inst.get<int64_t>(j));
            }
            qes::Instruction _inst(inst.get_name(), operands);
            // Copy annotations and properties over.
            for (annotation_t a : inst.get_annotations()) {
                _inst.put(a);
            }
            for (const auto& [ k, v ] : inst.get_property_map()) {
                _inst.put(k, v);
            }
            push_back_instruction(_inst);
        }
    } else {
        // Add the instruction to the graph. Add an edge to each ancestor.
        sptr<vertex_t> v = make_and_add_vertex(inst_ctr++);
        v->instruction = inst;
        if (inst.get_name() == "measure") {
            v->meas_order = meas_inst_ctr++;
        }
        for (size_t i = 0; i < k; i++) {
            int64_t q = inst.get<int64_t>(i);
            sptr<vertex_t> anc = get_ancestor_of(q);
            if (!contains(v, anc)) {
                make_and_add_edge(v, anc, false);
            }
            // Update ancestry:
            recent_ancestor_map[q] = v;
        }
    }
}

}   // graph
}   // qontra
