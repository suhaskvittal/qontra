/*
 *  author: Suhas Vittal
 *  date:   30 May 2024
 * */

#include "placc/fpn.h"
#include "placc/cx.h"

#include <qontra/tables.h>

#include <vtils/utility.h>

using namespace qontra;
using namespace graph;

static tables::ErrorAndTiming et;

namespace placc {

qes::Program<>
FPN::phase_one_schedule(fp_t& lat) {
    if (timestep_map.empty()) compute_cnot_order();
    qes::Program<> sch;
    CXManager cxm;

    std::cout << "[ p1 ] max timestep = " << max_timestep << std::endl;
    for (int mx = 0; mx <= 1; mx++) {
        lat += 2*et.t_g1q;
        push_back_gate(sch, "reset", {inplace_parity_qubits, flag_qubits});
        if (mx) push_back_gate(sch, "h", inplace_parity_qubits);
        else    push_back_gate(sch, "h", flag_qubits);
        // Entangle the flag qubits with the parity qubits according to the cnot_order.
        for (size_t t = 0; t <= max_timestep; t++) {
            for (sptr<fpn_v_t> f : flag_qubits) {
                if (!timestep_map[f].count(t)) continue;
                sptr<fpn_v_t> p = timestep_map[f][t];
                if (mx) cxm.push_back_cx(p->id, f->id);
                else    cxm.push_back_cx(f->id, p->id);
            }
        }
        lat += 2*cxm.get_depth()*et.t_g2q;
        auto cx_flag_ent = cxm.flush();
        vtils::push_back_range(sch, cx_flag_ent);
        // Perform data CX with flags.
        for (sptr<fpn_v_t> f : flag_qubits) {
            for (sptr<fpn_v_t> q : get_neighbors(f)) {
                if (q->qubit_type == fpn_v_t::type::data && !q->is_widowed) {
                    if (mx) cxm.push_back_cx(f->id,q->id);
                    else    cxm.push_back_cx(q->id,f->id);
                }
            }
        }
        lat += cxm.get_depth()*et.t_g2q;
        cxm.flush(sch);
        vtils::push_back_range(sch, cx_flag_ent);
        // Measure parity and flag qubits.
        lat += et.t_g1q;
        if (mx) push_back_gate(sch, "h", inplace_parity_qubits);
        else    push_back_gate(sch, "h", flag_qubits);
        lat += et.t_ro;
        push_back_gate(sch, "measure", {inplace_parity_qubits, flag_qubits});
    }
    return sch;
}

}   // placc
