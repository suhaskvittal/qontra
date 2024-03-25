/* author: Suhas Vittal
 *  date:   27 December 2023
 * */

#ifndef PROTEAN_NETWORK_h
#define PROTEAN_NETWORK_h

#include <qontra/ext/qes.h>
#include <qontra/ext/stim.h>
#include <qontra/graph/tanner_graph.h>

#include <vtils/bijective_map.h>
#include <vtils/two_level_map.h>

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

    bool is_check(void);
};

struct raw_edge_t : graph::base::edge_t {};

struct phys_vertex_t : graph::base::vertex_t {
    std::set<sptr<raw_vertex_t>> role_set;
    vtils::BijectiveMap<size_t, sptr<raw_vertex_t>> cycle_role_map;

    void consume(sptr<phys_vertex_t>);

    void add_role(sptr<raw_vertex_t>, size_t cycle);
    void delete_role(sptr<raw_vertex_t>);
    void push_back_role(sptr<raw_vertex_t>, size_t min_cycle=0);

    bool has_role_of_type(raw_vertex_t::type) const;
    
    void clear_roles(void);
private:
    std::map<raw_vertex_t::type, size_t> role_type_map;
};

struct phys_edge_t : graph::base::edge_t {
    size_t tsv_layer;

    bool is_out_of_plane(void) const;
};

}   // net

// The RawNetwork keeps track of the "roles", or all the qubit interactions that need to
// happen to implement syndrome extraction.
class RawNetwork : public graph::Graph<net::raw_vertex_t, net::raw_edge_t> {
public:
    struct parity_support_t {
        sptr<net::raw_vertex_t> check;

        std::set<sptr<net::raw_vertex_t>>   data;
        std::set<sptr<net::raw_vertex_t>>   flags;
        std::set<sptr<net::raw_vertex_t>>   proxies;

        std::set<sptr<net::raw_vertex_t>>   all;
    };

    RawNetwork(graph::TannerGraph*);
    RawNetwork(RawNetwork&&) = default;

    void delete_vertex(sptr<net::raw_vertex_t>) override;

    sptr<net::raw_vertex_t> make_and_add_vertex(void);

    // Input: data qubit, data qubit, parity qubit
    // Returns: flag qubit.
    sptr<net::raw_vertex_t> add_flag(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);
    // Input: xyz qubit, xyz qubit (xyz meaning can be any type).
    // Returns: proxy qubit. The proxy qubit splits the edge between the two input qubits.
    sptr<net::raw_vertex_t> add_proxy(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);
    // This does the same as the above command, but just takes an edge instead of two qubits.
    sptr<net::raw_vertex_t> add_proxy(sptr<net::raw_edge_t>);

    // This function merges the two roles and deletes one of the operands.
    // The merge operation moves the neighbors of the less-important role into 
    // the more-important role. The role hierarchy is as follows:
    //  proxy < flag < parity   (if one of the operands is a data qubit, there is undefined behavior).
    // If both roles are equally important, the second is merged into the first.
    //
    // Returns: the deleted vertex, or nullptr if nothing was deleted.
    sptr<net::raw_vertex_t> merge(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);
    
    // Input: two qubits (any type)
    // Output: a path (inclusive) of proxies between the two qubits. This is memoized.
    std::vector<sptr<net::raw_vertex_t>>&
        get_proxy_walk_path(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);
    // Input: a parity qubit
    // Output: the support of the parity qubit.
    parity_support_t&       get_support(sptr<net::raw_vertex_t>);
    // Input: two qubits (any type)
    // Output: returns the check they belong to, or returns nullptr if no check exists. This is memoized.
    sptr<net::raw_vertex_t> are_in_same_support(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);
    // This function is the same as above but extended for multiple qubits.
    // This is NOT memoized as this is a more niche use case.
    sptr<net::raw_vertex_t> are_in_same_support(std::initializer_list<sptr<net::raw_vertex_t>>);

    void    reset_memoization(void);
    void    disable_memoization(void);
    
    // Tanner graph tracking structures:
    vtils::BijectiveMap<sptr<graph::tanner::vertex_t>, sptr<net::raw_vertex_t>> v_tanner_raw_map;

    // Flag tracking structures:
    //
    // parity qubit --> list of flag qubits used in syndrome extraction of check.
    std::map<sptr<net::raw_vertex_t>, std::vector<sptr<net::raw_vertex_t>>> 
        flag_ownership_map;
    // parity qubit --> data qubit --> flag qubit used in syndrome extraction of check.
    vtils::TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>>
        flag_assignment_map;
    std::set<sptr<net::raw_vertex_t>> x_flag_set;

    // Scheduling structures:
    // 
    // parity qubit --> list of data qubits (or nullptr) in the order of syndrome extraction CNOTs.
    std::map<sptr<net::raw_vertex_t>, std::vector<sptr<net::raw_vertex_t>>>
        schedule_order_map;

    graph::TannerGraph* tanner_graph;
    bool enable_memoization; // This should be set to true once the network has stabilized.
private:
    void update_endpoint_in_indirection_map(
            sptr<net::raw_vertex_t> proxy,
            sptr<net::raw_vertex_t> old_endpoint,
            sptr<net::raw_vertex_t> new_endpoint);

    void flag_proxy_merge(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>);

    // This function replaces the first operand with the second in all tracking structures. A third operand,
    // where, allows the user to identify which tracking structures to replace in. where is a bitvector:
    //      bit 0 = flag_ownership_map
    //      bit 1 = flag_assignment_map
    void replace_in_tracking_structures(sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t> with, uint8_t where=7);

    vtils::TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, std::vector<sptr<net::raw_vertex_t>>>
        proxy_memo_map;
    std::map<sptr<net::raw_vertex_t>, parity_support_t>
        support_memo_map;
    vtils::TwoLevelMap<sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>, sptr<net::raw_vertex_t>>
        same_support_memo_map;

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
    PhysicalNetwork(graph::TannerGraph*);
    PhysicalNetwork(PhysicalNetwork&&) = default;
    
    // All of the graph modification functions need to be updated to handle planarity.
    // Note that deletes only free up spots on the processor bulk. Only adding edges
    // requires checking for planarity.
    bool add_vertex(sptr<net::phys_vertex_t>) override;
    bool add_edge(sptr<net::phys_edge_t>) override;
    void delete_vertex(sptr<net::phys_vertex_t>) override;
    void delete_edge(sptr<net::phys_edge_t>) override;

    sptr<net::phys_vertex_t> make_and_add_vertex(void);

    // If the edge can be moved to the immediately lower processor layer, it is done and the
    // function returns true. Otherwise, it returns false.
    bool test_and_move_edge_down(sptr<net::phys_edge_t>); 

    void load_colors_from_file(std::string);
    // Assigns colors to the checks greedily. This is used for memory experiments of color
    // codes.
    void assign_colors_to_checks(void);

    size_t  get_thickness(void);
    size_t  get_bulk_degree(sptr<net::phys_vertex_t>);
    fp_t    get_round_latency(void);
    size_t  get_round_cnots(void);
    fp_t    get_expected_collisions(void);

    graph::TannerGraph* get_tanner_graph(void);
    uptr<RawNetwork>&   get_raw_connection_network(void);

    sptr<net::phys_vertex_t> get_physical_qubit_for(sptr<net::raw_vertex_t>);

    // There are the passes for Protean. Each pass returns true if some change was made.
    //
    // Combines physical qubits with identical connectivity (i.e. seen in color codes).
    bool join_qubits_with_identical_support(void);
    // Same as above, but when the neighbors of one are a strict subset of another.
    bool join_qubits_with_partial_support(void);
    // Adds flag qubits to the architecture.
    bool make_flags(void);
    // Reduces connectivity by adding proxy qubits.
    bool add_connectivity_reducing_proxies(void);
    // Removes qubits with degree-1 or degree-2 connectivity.
    bool contract_small_degree_qubits(void);
    // Try and place out-of-plane edges into the processor bulk, and deletes any
    // empty processor layers.
    bool reallocate_edges(void);
    // Recomputes the cycle_role_map of each phys_vertex_t. During the optimizations, the cycle_role_map
    // is a partial ordering of roles. However, to get the total order, we must consider edges between
    // roles as well. This function refreshes the maps to reflect this ordering.
    bool recompute_cycle_role_maps(void);
    // Relabel qubits so that the ids are not sporadic.
    bool relabel_qubits(void);

    // Runs the scheduler and yield simulator and updates any variables/statistics.
    void finalize(void);

    // IO:
    void    write_stats_file(std::string);
    void    write_schedule_folder(std::string);
    void    write_coupling_file(std::string);
    void    write_role_file(std::string);
    void    write_tanner_graph_file(std::string);
    void    write_flag_assignment_file(std::string);

    struct {
        size_t max_connectivity = 4;
        size_t max_thickness = 1;   // 0 = means only processor bulk, n = n TSV layers.

        fp_t min_qubit_frequency = 5.00e9;
        fp_t max_qubit_frequency = 5.34e9;
        fp_t qubit_frequency_step = 0.01e9;
        fp_t fabrication_precision = 30e6;

        // Scheduler settings.
        size_t  rounds = 1;

        // Optimization settings.
        bool    force_unopt_flags = false;
        bool    force_xz_flag_merge = false;
        bool    enable_flag_reduction = false;
        bool    enable_proxy_triangle_search = false;
    } config;
private:
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

    // Allocates a new processor layer.
    uptr<ProcessorLayer>& push_back_new_processor_layer(void);
    // Consumes a physical qubit safely.
    void consume(sptr<net::phys_vertex_t>, sptr<net::phys_vertex_t>);
    // Returns true if two physical qubits are contained within the same support for any pair of roles.
    bool are_in_same_support(sptr<net::phys_vertex_t>, sptr<net::phys_vertex_t>);

    // The reason tanner_graph is a raw pointer whereas raw_connection_network is a unique ptr
    // is that tanner_graph is supplied by the user. The user may choose to allocate it on the stack.
    graph::TannerGraph* tanner_graph;

    // raw_connection_network contains all roles in the network, from proxy to flag to data (etc.).
    // Each phys_vertex_t corresponds to at least one raw_vertex_t (if not more).
    uptr<RawNetwork> raw_connection_network;
    std::map<sptr<net::raw_vertex_t>, sptr<net::phys_vertex_t>> role_to_phys;

    // processor_layers contains the physical placement of edges in the processor. processor_layers[0]
    // always corresponds to the processor bulk (lowest layer), and other layers are the TSV layers.
    std::vector<uptr<ProcessorLayer>> processor_layers;

    // Other tracking structures:
    std::map<sptr<net::raw_vertex_t>, int> check_color_map;
    // Id tracking for make_and_add_vertex 
    uint64_t id_ctr;

    // Other statistics:
    fp_t round_latency;
    size_t round_cnots;
    fp_t expected_collisions;

    qes::Program<> x_memory;
    qes::Program<> z_memory;

    friend class Scheduler;
};

}   // protean
}   // qontra

#include "network.inl"

#endif  // PROTEAN_NETWORK_h
