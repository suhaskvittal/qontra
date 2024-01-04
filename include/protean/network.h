/* 
 *  author: Suhas Vittal
 *  date:   27 December 2023
 * */

#ifndef PROTEAN_NETWORK_h
#define PROTEAN_NETWORK_h

#include <defs/bijective_map.h>
#include <defs/two_level_map.h>
#include <graph/tanner_graph.h>
#include <instruction.h>
#include <stimext.h>

namespace qontra {
namespace protean {

// We will have two types of networks:
//  (1) A "raw" network, which has the following rules:
//      (i) a qubit can only be used for "one" role (i.e. parity qubit for a single check, a flag, etc.)
//  (2) A physical network, which is the actual coupling graph. Here, qubits can have multiple roles
//      at different points in time.

namespace net {

struct raw_vertex_t : graph::base::vertex_t {
    enum class type { data, xparity, zparity, flag, proxy };
    type qubit_type;
};

struct raw_edge_t : graph::base::edge_t {};

struct phys_vertex_t : graph::base::vertex_t {
    std::set<sptr<raw_vertex_t>> role_set;
    BijectiveMap<size_t, sptr<raw_vertex_t>> cycle_role_map;

    void consume(sptr<phys_vertex_t>);

    void add_role(sptr<raw_vertex_t>, size_t cycle);
    void push_back_role(sptr<raw_vertex_t>, size_t min_cycle=0);
    bool has_role_of_type(raw_vertex_t::type);
private:
    std::set<raw_vertex_t::type> role_type_set;
};

struct phys_edge_t : graph::base::edge_t {
    size_t tsv_layer;

    bool is_out_of_plane(void);
};

}   // net

// The RawNetwork keeps track of the "roles", or all the qubit interactions that need to
// happen to implement syndrome extraction.
class RawNetwork : public graph::Graph<net::raw_vertex_t, net::raw_edge_t> {
public:
    RawNetwork(graph::TannerGraph&);

    sptr<net::raw_vertex_t> make_vertex(void) override;

    // Input: data qubit, data qubit, parity qubit
    // Returns: flag qubit.
    sptr<net::raw_vertex_t> add_flag(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);
    // Input: xyz qubit, xyz qubit (xyz meaning can be any type).
    // Returns: proxy qubit. The proxy qubit splits the edge between the two input qubits.
    sptr<net::raw_vertex_t> add_proxy(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);
    // This does the same as the above command, but just takes an edge instead of two qubits.
    sptr<net::raw_vertex_t> add_proxy(sptr<net::raw_edge_t>);

    // Input: xyz qubit, proxy qubit, a reference to a vector which will be populated with the walk.
    // Returns: non-proxy qubit obtained by walking through the proxy_indirection_map.
    sptr<net::raw_vertex_t> proxy_walk(sptr<net::raw_vertex_t>,
                                            sptr<net::raw_vertex_t>,
                                            std::vector<sptr<net::raw_vertex_t>>&);
    
    // Tanner graph tracking structures:
    BijectiveMap<sptr<graph::tanner::vertex_t>, sptr<net::raw_vertex_t>> v_tanner_raw_map;

    // Flag and proxy tracking structures:
    //
    // parity qubit --> list of flag qubits used in syndrome extraction of check.
    std::map<sptr<net::raw_vertex_t>, std::vector<sptr<net::raw_vertex_t>>> 
        flag_ownership_map;
    // parity qubit --> data qubit --> flag qubit used in syndrome extraction of check.
    TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>>
        flag_assignment_map;
    // proxy qubit --> xyz qubit --> qubit xyz was linked to before
    TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>>
        proxy_indirection_map;

    // Scheduling structures:
    // 
    // parity qubit --> list of data qubits (or nullptr) in the order of syndrome extraction CNOTs.
    std::map<sptr<net::raw_vertex_t>, std::vector<sptr<net::raw_vertex_t>>>
        schedule_order_map;

    graph::TannerGraph tanner_graph;
private:
    uint64_t id_ctr = 0;
};

// The ProcessorLayer class encapsulates a layer of the processor that is planar.
// An array of these layers builds a full processor.
class ProcessorLayer : public graph::Graph<net::phys_vertex_t, net::phys_edge_t> {
public:
    bool add_edge(sptr<net::phys_edge_t>) override;

    size_t get_max_endpoint_degree(sptr<net::phys_edge_t>);

    bool is_planar(void);
protected:
    bool update_state(void) override;
private:
    bool _is_planar; // Never access this variable directly!!!
};

// A PhysicalNetwork is the realization of the processor.
class PhysicalNetwork : public graph::Graph<net::phys_vertex_t, net::phys_edge_t> {
public:
    PhysicalNetwork(graph::TannerGraph&);
    
    static PhysicalNetwork from_folder(std::string);

    // All of the graph modification functions need to be updated to handle planarity.
    // Note that deletes only free up spots on the processor bulk. Only adding edges
    // requires checking for planarity.
    bool add_vertex(sptr<net::phys_vertex_t>) override;
    bool add_edge(sptr<net::phys_edge_t>) override;
    void delete_vertex(sptr<net::phys_vertex_t>) override;
    void delete_edge(sptr<net::phys_edge_t>) override;

    sptr<net::phys_vertex_t> make_vertex(void) override;

    // If the edge can be moved to the immediately lower processor layer, it is done and the
    // function returns true. Otherwise, it returns false.
    bool test_and_move_edge_down(sptr<net::phys_edge_t>); 

    size_t  get_thickness(void);
    size_t  get_bulk_degree(sptr<net::phys_vertex_t>);

    // Optimizations:
    //
    //  join_qubits_with_identical_support
    //      -- merges qubits connected to the same data qubits (i.e. as in the color code).
    //  join_qubits_with_partial_support
    //      -- similar to above, but requires that the adjacency list is a strict subset.
    //  make_flags
    //      -- creates flag qubits and fault-tolerant circuits.
    //      -- PRECONDITION: data qubits are still connected directly to parity qubits.
    //  add_connectivity_reducing_proxies
    //      -- creates proxy qubits wherever the connectivity of a qubit is too high.
    //      -- proxies are not guaranteed to be fault-tolerant when connected to data qubits,
    //          so only use them after adding flag qubits.
    //  contract_small_degree_qubits
    //      -- contracts qubits (that are non-data!!!) that are connected to at most two neighbors,
    //          as this contraction does not create new connectivity violations.
    //  reallocate_edges
    //      -- attempts to move edges down processor layers, and deletes any empty layers
    //  relabel_qubits
    //      -- this is more of a convenience pass
    //      -- It may be that we can have missing ids during execution due to modifications. This
    //          pass fixes that. Recommended to run at the END of compilation.
    //
    //  All passes return true if some modification was made,
    bool join_qubits_with_identical_support(void);
    bool join_qubits_with_partial_support(void);
    bool make_flags(void);
    bool add_connectivity_reducing_proxies(void);
    bool contract_small_degree_qubits(void);
    // Try and place out-of-plane edges into the processor bulk, and deletes any
    // empty processor layers.
    bool reallocate_edges(void);
    // Relabel qubits so that the ids are not sporadic.
    bool relabel_qubits(void);

    RawNetwork get_raw_connection_network(void);

    schedule_t make_schedule(void);

    struct {
        size_t max_connectivity = 4;
        size_t max_thickness = 1;   // 0 = means only processor bulk, n = n TSV layers.
    } config;
private:

    // Allocates a new processor layer.
    ProcessorLayer& push_back_new_processor_layer(void);

    void consume(sptr<net::phys_vertex_t>, sptr<net::phys_vertex_t>);

    // Determines if a proposed flag (represented by the two raw vertices -- these are data qubits)
    // is actually useful -- that it protects against a weight-2 error.
    //
    // To do so, we will use an operator tree, which computes logical operators via dynamic programming.
    // 
    // As the traversal is at worst exponential time in the size of the code, we will just do all flags
    // at once.
    stim::simd_bits<SIMD_WIDTH> do_flags_protect_weight_two_error(
            std::set<std::pair<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>>>, 
            bool is_x_error);

    // raw_connection_network contains all roles in the network, from proxy to flag to data (etc.).
    // Each phys_vertex_t corresponds to at least one raw_vertex_t (if not more).
    RawNetwork raw_connection_network;
    std::map<sptr<net::raw_vertex_t>, sptr<net::phys_vertex_t>> role_to_phys;
    // processor_layers contains the physical placement of edges in the processor. processor_layers[0]
    // always corresponds to the processor bulk (lowest layer), and other layers are the TSV layers.
    std::vector<ProcessorLayer> processor_layers;
    // Other tracking structures:
    //
    // Tracks the heights of TSV edges for each vertex.
    std::map<sptr<net::phys_vertex_t>, std::set<size_t>> occupied_tsvs;
    
    uint64_t id_ctr;

    friend void write_schedule_file(std::string, PhysicalNetwork&);
    friend void write_coupling_file(std::string, PhysicalNetwork&);
    friend void write_role_file(std::string, PhysicalNetwork&);
    friend void write_tanner_graph_file(std::string, PhysicalNetwork&);
    friend void write_flag_assignment_file(std::string, PhysicalNetwork&);
};

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

void write_schedule_file(std::string, PhysicalNetwork&);
void write_coupling_file(std::string, PhysicalNetwork&);
void write_role_file(std::string, PhysicalNetwork&);
void write_tanner_graph_file(std::string, PhysicalNetwork&);
void write_flag_assignment_file(std::string, PhysicalNetwork&);

}   // protean
}   // qontra

#include "network.inl"

#endif  // PROTEAN_NETWORK_h
