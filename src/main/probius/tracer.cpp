/*
 *  author: Suhas Vittal
 *  date:   10 December 2023
 * */

#include <parsing/cmd.h>
#include <sim/memory_sim.h>
#include <tables.h>

#include <mpi.h>

using namespace qontra;

const std::map<std::string, MemorySimulator::lrc_policy_t>
LRC_POLICY_TABLE {
    std::make_pair("always", MemorySimulator::lrc_policy_t::always),
    std::make_pair("optimal", MemorySimulator::lrc_policy_t::optimal),
    std::make_pair("eraser", MemorySimulator::lrc_policy_t::eraser)
};

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    CmdParser pp(argc, argv);

    uint32_t distance;
    uint32_t rounds;
    uint64_t shots;
    fp_t error_rate;

    std::string stats_output_file;
    std::string stim_output_file;
    std::string syndrome_output_folder;

    bool is_memory_x = false;
    bool sim_leakage = false;

    std::string lrc_policy_name;
    MemorySimulator::lrc_policy_t lrc_policy = MemorySimulator::lrc_policy_t::none;
    MemorySimulator::lrc_circuit_t lrc_circuit = MemorySimulator::lrc_circuit_t::swap;

    //
    // Get command line options
    //

    // Basic options:

    if (!pp.get_uint32("d", distance)) return 1;
    if (!pp.get_float("p", error_rate)) return 1;
    if (!pp.get_uint64("shots", shots)) return 1;

    rounds = distance;
    pp.get_uint32("rounds", rounds);

    if (!pp.get_string("stats-out", stats_output_file)) return 1;
    if (!pp.get_string("stim-out", stim_output_file)) return 1;
    if (!pp.get_string("trace-out", syndrome_output_folder)) return 1;
    
    is_memory_x = pp.option_set("mx");

    // Leakage options:

    sim_leakage = pp.option_set("enl");
    if (pp.get_string("lrc-policy", lrc_policy_name)) lrc_policy = LRC_POLICY_TABLE.at(lrc_policy_name);
    if (pp.option_set("dqlr")) lrc_circuit = MemorySimulator::lrc_circuit_t::dqlr;

    //
    // Instantiate simulator and setup config.
    //

    configure_optimal_batch_size();
    if (world_rank == 0) {
        std::cout << "\tBatch size is " << experiments::G_SHOTS_PER_BATCH << "\n";
    }
    
    graph::LatticeGraph lattice = sim::surface_code_lattice_graph(distance);
    MemorySimulator sim(lattice);
    sim.config.lrc_policy = lrc_policy;
    sim.config.lrc_circuit = lrc_circuit;

    // Instantiate error model.
    const uint n = lattice.get_vertices().size();
    tables::ErrorAndTiming et;

    et.e_g1q = 1e-4;
    et.e_g2q = 1e-3;
    et.e_m1w0 = 1e-3;
    et.e_m0w1 = 3e-3;
    et.e_idle = 1e-4;

    et.t1 = 500e3;
    et.t2 = 250e3;

    et = et * (1000*error_rate);
    tables::populate(n, sim.config.errors, sim.config.timing, et);
    // Add leakage errors.
    if (sim_leakage) {
        sim.config.enable_leakage = true;
        for (uint i = 0; i < n; i++) {
            for (uint j = 0; j < n; j++) {
                if (i == j) continue;
                auto i_j = std::make_pair(i, j);
                sim.config.errors.op2q_leakage_transport["cx"][i_j] = 0.1;
                sim.config.errors.op2q_leakage_injection["cx"][i_j] = 0.1*error_rate;
                sim.config.errors.op2q_leakage_transport["liswap"][i_j] = 0.1;
                sim.config.errors.op2q_leakage_injection["liswap"][i_j] = 0.1*error_rate;
            }
        }
    }
    // Set reset error to p.
    for (uint i = 0; i < n; i++) {
        sim.config.errors.op1q["reset"][i] = error_rate;
    }

    sim.config.distance = distance;
    sim.config.rounds = rounds;
    sim.config.stim_output_file = stim_output_file;
    sim.config.syndrome_output_folder = syndrome_output_folder;
    sim.config.data_output_file = stats_output_file;
    sim.config.is_memory_x = is_memory_x;

    sim.run(shots);

    MPI_Finalize();
    return 0;
}
