/*
 *  author: Suhas Vittal
 *  date:   4 March 2024
 * */

#include <qontra/experiments/memory.h>
#include <qontra/ext/stim.h>
#include <qontra/graph/decoding_graph.h>

#include <stdlib.h>

using namespace qontra;
using namespace graph;
using namespace decoding;

std::set<sptr<hyperedge_t>>
collect_degenerate_errors(uptr<DecodingGraph>& gr) {
    std::set<sptr<hyperedge_t>> degens;
    for (sptr<hyperedge_t> e : gr->get_all_edges()) {
        if (e->get_order() == 0) {
            degens.insert(e);
            continue;
        }
        // Otherwise, check the equivalence class of the edge.
        EdgeClass ec = gr->get_edge_class(e);
        sptr<hyperedge_t> r = ec.get_representative();
        if (r->flags == e->flags && r->frames != e->frames) {
            degens.insert(r);
            degens.insert(e);
        }
    }
    return degens;
}

int main(int argc, char* argv[]) {
    std::string qes_file(argv[1]);
    int flips_per_error = atoi(argv[2]);

    DetailedStimCircuit circuit = make_circuit(qes_file, 1e-3);
    uptr<DecodingGraph> gr = std::make_unique<DecodingGraph>(circuit, flips_per_error);
    std::set<sptr<hyperedge_t>> degens = collect_degenerate_errors(gr);

    std::cout << circuit << std::endl;

    // Degenerate detector errors.
    stim::DetectorErrorModel ddem;
    for (sptr<hyperedge_t> e : degens) {
        std::vector<stim::DemTarget> targets;
        for (sptr<vertex_t> v : e->get<vertex_t>()) {
            if (v->is_boundary_vertex) continue;
            targets.push_back(stim::DemTarget::relative_detector_id(v->id));
        }
        for (uint64_t f : e->flags) {
            targets.push_back(stim::DemTarget::relative_detector_id(f));
        }
        for (uint64_t fr : e->frames) {
            targets.push_back(stim::DemTarget::observable_id(fr));
        }
        stim::SpanRef<stim::DemTarget> sp(&targets[0], &targets[0] + targets.size());
        ddem.append_error_instruction(e->probability, sp);
    }
    auto explained_errors = stim::ErrorMatcher::explain_errors_from_circuit(circuit, &ddem, false);
    for (auto ee : explained_errors) {
        std::cout << ee << std::endl;
    }
    return 0;
}
