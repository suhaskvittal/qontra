/*
 * author:  Suhas Vittal
 * date:    23 May 2024
 * */

#include "gen.h"

#include "qontra/decoder.h"
#include "qontra/decoder/concat_mwpm.h"
#include "qontra/decoder/restriction.h"
#include "qontra/decoder/mobius.h"
#include "qontra/graph/decoding_graph.h"

#include "qontra/experiments/memory.h"

#include <random>

#include <stdlib.h>

using namespace qontra;
using namespace graph;
using namespace decoding;

using namespace vtils;

static std::mt19937_64 rng(0);

void
sample_from_graph(
        DecodingGraph& gr,
        const std::vector<sptr<hyperedge_t>> edges,
        int error_weight,
        stim::simd_bits_range_ref<SIMD_WIDTH> syndrome,
        stim::simd_bits_range_ref<SIMD_WIDTH> obs) 
{
    const size_t m = edges.size();

    sptr<hyperedge_t> e = edges.at(rng() % m);
    std::set<sptr<hyperedge_t>> visited;
    while (error_weight--) {
        for (sptr<vertex_t> v : e->get<vertex_t>()) {
            if (v->is_boundary_vertex) {
                continue;
            }
            syndrome[v->id] ^= 1;
        }
        for (uint64_t f : e->flags)                 syndrome[f] = 1;
        for (uint64_t fr : e->frames)               obs[fr] ^= 1;
        EdgeClass cl = gr.get_edge_class(e);
        visited.insert(cl.get_representative());
        
        std::cout << "Injected error [";
        for (sptr<vertex_t> v : e->get<vertex_t>()) std::cout << " " << print_v(v);
        std::cout << " ], flags = [";
        for (uint64_t f : e->flags) std::cout << " " << f;
        std::cout << " ], frames =";
        for (uint64_t fr : e->frames) std::cout << " " << fr;
        std::cout << std::endl;

        if (error_weight == 0) break;
        // Find the next candidate edge.
        std::vector<sptr<hyperedge_t>> candidates;
        for (sptr<hyperedge_t> x : edges) {
            if (x->flags.size()) continue;
            EdgeClass _cl = gr.get_edge_class(x);
            if (visited.count(_cl.get_representative())) continue;
            size_t common = 0;
            for (sptr<vertex_t> v : e->get<vertex_t>()) {
                if (std::find(x->endpoints.begin(), x->endpoints.end(), v) != x->endpoints.end()) {
                    common++;
                }
            }
            if (common > 0) candidates.push_back(x);
        }
        e = candidates[rng() % candidates.size()];
    }
}

int main(int argc, char* argv[]) {
    std::string tanner_graph_file(argv[1]);
    int error_weight = atoi(argv[2]);
    uint64_t shots = atoll(argv[3]);

    TannerGraph tgr = 
        create_graph_from_file<TannerGraph>(tanner_graph_file, io::update_tanner_graph);
    std::map<sptr<tanner::vertex_t>, int> color_map;
    if (tgr.compute_check_color_map(color_map) > 2) {
        std::cerr << "Found >3-coloring." << std::endl;
        exit(1);
    }
    DetailedStimCircuit circuit = make_capacity(&tgr, 1e-3, false, color_map);
    uptr<Decoder> dec = std::make_unique<ConcatMWPMDecoder>(circuit);
    // The number of errors does not matter as we don't care about boundaries.
    DecodingGraph gr(circuit, 1000);
    auto edges = gr.get_all_edges();
    // Remove all order 0 edges from the list.
    for (auto it = edges.begin(); it != edges.end(); ) {
        if ((*it)->get_order() == 0) it = edges.erase(it);
        else it++;
    }
    const size_t n_det = circuit.count_detectors(),
                 n_obs = circuit.count_observables();
    for (uint64_t s = 0; s < shots; s++) {
        std::cout << s << "--------------------------------------" << std::endl;
        stim::simd_bits<SIMD_WIDTH> syndrome(n_det),
                                    obs(n_obs);
        sample_from_graph(gr, edges, error_weight, syndrome, obs);
        auto res = dec->decode_error(syndrome);
        if (res.corr != obs) {
            std::cout << "is logical error!" << std::endl;
            std::cout << "\texpected: ";
            for (size_t i = 0; i < n_obs; i++) std::cout << obs[i]+0;
            std::cout << std::endl << "\treceived: ";
            for (size_t i = 0; i < n_obs; i++) std::cout << res.corr[i]+0;
            std::cout << std::endl;
        }
    }
    return 0;
}
