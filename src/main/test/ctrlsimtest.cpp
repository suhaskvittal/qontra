/*
 *  author: Suhas Vittal
 *  date:   23 June 2023
 * */

#include "decoder/mwpm.h"
#include "instruction.h"
#include "parsing/cmd.h"
#include "sim/control_sim.h"

#include <stim.h>

#include <iostream>
#include <vector>

using namespace qontra;

int main(int argc, char* argv[]) {
    CmdParser cmd_parser(argc, argv);

    schedule_t  prog;
    uint        n;

    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (cmd_parser.option_set("asm")) {
        std::string asm_file;
        cmd_parser.get_string("asm", asm_file);
        prog = from_file(asm_file);

        cmd_parser.get_uint32("size", n);
    } else {
        stim::CircuitGenParameters circ_params(3, 3, "rotated_memory_z");
        auto circ = stim::generate_surface_code_circuit(circ_params).circuit;

        n = from_stim_circuit(circ, prog);

        prog.push_back({"decode", {0}, {}});
        prog.push_back({"xorfr", {0, 0}, {}});
        prog.push_back({"savem", {}, {}});
        prog.push_back({"done", {}, {}});
    }

    experiments::G_SHOTS_PER_BATCH = 1 << 12;

    ControlSimulator sim(n, prog);
    tables::ErrorAndTiming params;
    tables::populate(n, sim.params.errors, sim.params.timing, params);
    sim.build_canonical_circuit();

    stim::Circuit error_model = sim.get_canonical_circuit();
    decoder::MWPMDecoder mwpm(error_model);

    sim.load_decoder(&mwpm);
    
    // Check expected LER.
    auto mxp_res = memory_experiment(&mwpm, 1L << 20);
    if (world_rank == 0) {
        std::cout << "Expected LER = " << mxp_res.logical_error_rate << "\n";
    }

    sim.params.verbose = cmd_parser.option_set("v") && (world_rank == 0);

    sim.run(1L << 20);
    if (world_rank == 0) {
        std::cout << "Probability histogram:\n";
        for (auto pair : sim.prob_histograms) {
            std::cout << "\t" << std::hex;
            for (uint64_t x : pair.first)   std::cout << x << " ";
            std::cout << "\t" << std::dec << pair.second << "\n";
        }
        std::cout << "Time taken to simulate: " << sim.sim_time * 1e-9 << "\n";
    }
    MPI_Finalize();
}

