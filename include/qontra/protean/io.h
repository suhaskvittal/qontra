/*
 *  author: Suhas Vittal
 *  date:   12 January 2024
 * */

#ifndef PROTEAN_IO_h
#define PROTEAN_IO_h

#include "qontra/protean/network.h"

namespace qontra {
namespace protean {

// Takes in a string of characters corresponding to pass execution and updates the PhysicalNetwork
// reference accordingly.
//
// Each pass is a three letter character:
//      Jid -- join_qubits_with_identical_support
//      Jpa -- join_qubits_with_partial_support
//      Fla -- make_flags
//      Prx -- add_connectivity_reducing_proxies
//      Con -- contract_small_degree_qubits
//      Ral -- reallocate_qubits
//      Rcr -- recompute_cycle_role_maps
//      Rlb -- relabel_qubits
// The above is not case sensitive.
//
// Special characters:
//      .,;:| -- separates passes out for readability (ignored)
//      (...) -- the passes in the parenthesis are a single group of passes
//      PASS+, (...)+ -- executes the PASS or pass group until no modifications are made
//  Example string:
//      Jid.Ral.Fla.Ral.Jpa.Ral.(Prx.Con.Jpa.Ral)+.Rlb.
//
//      Equivalents (examples):
//          JidRalFlaRalJpaRal(PrxConJpaRal)+Rlb
//          jidralflaraljparal(prxconjparal)+rlb
//          jid|ral|fla|ral|jpa|ral(prx.con.jpa.ral)+rlb
//      Use whichever is most readable.
bool update_network(std::string pass_string, PhysicalNetwork&, bool verbose=false);

// Writes the physical network to a folder.
// Generated Files:
//  (1) round.asm: ASM for a round of syndrome extraction
//  (2) coupling.txt: coupling graph for network
//  (2) roles.txt: corresponding roles and cycles for each (physical) parity qubit.
//  (3) tanner_graph.txt: corresponding checks for each (raw) parity qubit
//                      and logical observables (not directly usable with TannerGraph)
//  (4) flag_assignment.txt: corresponding flags for each (raw) parity qubit
void write_network_to_folder(std::string, PhysicalNetwork&);

void write_stats_file(std::string, PhysicalNetwork&);
void write_schedule_file(std::string, PhysicalNetwork&);
void write_coupling_file(std::string, PhysicalNetwork&);
void write_role_file(std::string, PhysicalNetwork&);
void write_tanner_graph_file(std::string, PhysicalNetwork&);
void write_flag_assignment_file(std::string, PhysicalNetwork&);

}   // protean
}   // qontra

#include "io.inl"

#endif  // PROTEAN_IO_h
