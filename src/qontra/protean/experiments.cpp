/*
 *  author: Suhas Vittal
 *  date:   2 March 2024
 * */

#include "qontra/protean/experiments.h"

#include <fstream>
#include <iostream>

#include <math.h>

namespace qontra {
namespace protean {

using namespace graph;
using namespace base;

inline fp_t get_cx_error(uptr<CouplingGraph>& graph, sptr<edge_t> e) {
    sptr<vertex_t> v = e->get_source<vertex_t>(),
                    w = e->get_target<vertex_t>();
    const fp_t dgs = 0.5*static_cast<fp_t>( graph->get_degree(v) + graph->get_degree(w) );
    // If dgs > 4.0, then e^(dgs - 4.0).
    // else 1 + 0.1*(dgs-4.0)
    return dgs >= 4.0 ? pow(M_E, dgs-4.0) : 1.0 + 0.15*(dgs-4.0);
}

void
safe_add_edge(uptr<CouplingGraph>& graph, uint64_t src, uint64_t dst) {
    sptr<vertex_t> v, w;
    if (graph->contains(src)) {
        v = graph->get_vertex(src);
    } else {
        v = graph->make_and_add_vertex(src);
    }
    if (graph->contains(dst)) {
        w = graph->get_vertex(dst);
    } else {
        w = graph->make_and_add_vertex(dst);
    }
    graph->make_and_add_edge(v, w);
}

uptr<CouplingGraph>
read_coupling_graph(std::string coupling_file) {
    uptr<CouplingGraph> gr = std::make_unique<CouplingGraph>();
    
    std::ifstream fin(coupling_file);
    std::string ln;

    // Dump the first line -- we don't care about it (just a header for readability).
    std::getline(fin, ln);
    while (std::getline(fin, ln)) {
        size_t pos1, pos2;
        pos1 = ln.find(",");
        pos2 = ln.find(",", pos1+1);

        std::string s1 = ln.substr(0, pos1),        // First vertex
                    s2 = ln.substr(pos1+1, pos2),   // Second vertex
                    s3 = ln.substr(pos2+1);         // Processor Layer (currently unused).
        uint64_t v = std::stoull(s1),
                 w = std::stoull(s2);
        safe_add_edge(gr, v, w);
    }
    return gr;
}

void
make_error_and_timing_from_coupling_graph(std::string coupling_file, ErrorTable& errors, TimeTable& timing) {
    uptr<CouplingGraph> gr = std::move(read_coupling_graph(coupling_file));

    for (sptr<vertex_t> v : gr->get_vertices()) {
        const uint64_t q = v->id;

        errors.op1q["h"][q] = 0.1;
        errors.op1q["reset"][q] = 0.1;
        errors.idling[q] = 0.0;//0.1;
        errors.m1w0[q] = 1.0;
        errors.m0w1[q] = 1.0;

        timing.op1q["h"][q] = 30;
        timing.op1q["reset"][q] = 30;
        timing.op1q["measure"][q] = 800;
        timing.t1[q] = 300;
        timing.t2[q] = 150;
    }

    for (sptr<edge_t> e : gr->get_edges()) {
        sptr<vertex_t> v = e->get_source<vertex_t>(),
                        w = e->get_target<vertex_t>();
        const uint64_t q1 = v->id,
                         q2 = w->id;
        fp_t cxp = 0.0;//get_cx_error(gr, e);
        // Populate table.
        auto q1_q2 = std::make_pair(q1, q2),
             q2_q1 = std::make_pair(q2, q1);
        errors.op2q["cx"][q1_q2] = cxp;
        errors.op2q["cx"][q2_q1] = cxp;
        timing.op2q["cx"][q1_q2] = 40;
        timing.op2q["cx"][q2_q1] = 40;
    }
}


}   // protean
}   // qontra
