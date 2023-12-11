/*
 *  author: Suhas Vittal
 *  date:   10 December 2023
 * */

#include <parsing/cmd.h>
#include <sim/memory_sim.h>
#include <tables.h>

#include <mpi.h>

using namespace qontra;

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);

    CmdParser pp(argc, argv);

    uint32_t distance;
    uint32_t rounds;
    uint64_t shots;
    fp_t error_rate;

    std::string stim_output_file;
    std::string syndrome_output_folder;

    bool is_memory_x = false;

    if (!pp.get_uint32("d", distance)) return 1;
    if (!pp.get_float("p", error_rate)) return 1;
    if (!pp.get_uint64("shots", shots)) return 1;

    rounds = distance;
    pp.get_uint32("rounds", rounds);

    if (!pp.get_string("stim-out", stim_output_file)) return 1;
    if (!pp.get_string("trace-out", syndrome_output_folder)) return 1;
    
    is_memory_x = pp.option_set("mx");

    // Make trace folder if it does not exist.
    std::filesystem::path syndrome_output_folder_path(syndrome_output_folder);
    if (!std::filesystem::exists(syndrome_output_folder_path)) {
        safe_create_directory(syndrome_output_folder_path);
    }

    graph::LatticeGraph lattice = sim::surface_code_lattice_graph(distance);
    MemorySimulator sim(lattice);

    // Instantiate error model.
    tables::ErrorAndTiming et;
    et = et * (1000*error_rate);
    tables::populate(lattice.get_vertices().size(), sim.config.errors, sim.config.timing, et);

    sim.config.distance = distance;
    sim.config.rounds = rounds;
    sim.config.stim_output_file = stim_output_file;
    sim.config.syndrome_output_folder = syndrome_output_folder;
    sim.config.is_memory_x = is_memory_x;

    sim.run(shots);

    MPI_Finalize();
    return 0;
}
