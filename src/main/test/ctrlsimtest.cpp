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
    CmdParser parser(argc, argv);

    schedule_t  prog;
    uint        n;

    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (parser.option_set("asm")) {
        std::string asm_file;
        parser.get_string("asm", asm_file);
        prog = schedule_from_file(asm_file);

        parser.get_uint32("size", n);
    } else {
        stim::CircuitGenParameters circ_params(3, 3, "rotated_memory_z");
        auto circ = stim::generate_surface_code_circuit(circ_params).circuit;

        prog = schedule_from_stim(circ);
        n = circ.count_qubits();

        prog.push_back({"decode", {0}, {}});
        prog.push_back({"xorfr", {0, 0}, {}});
        prog.push_back({"savem", {}, {}});
        prog.push_back({"done", {}, {}});
    }

    uint64_t shots;
    if (!parser.get_uint64("shots", shots)) return 1;

    experiments::G_SHOTS_PER_BATCH = 1 << 12;

    FrameSimulator fsim(n, experiments::G_SHOTS_PER_BATCH);

    ControlSimulator sim(n, prog, &fsim);
    tables::ErrorAndTiming params;
    tables::populate(n, sim.params.errors, sim.params.timing, params);

    decoder::MWPMDecoder* mwpm = nullptr;
    if (!parser.option_set("nodec")) {
        sim.build_error_model();
        stim::Circuit error_model = sim.get_error_model();
        mwpm = new decoder::MWPMDecoder(error_model);

        // Check expected LER.
        auto mxp_res = memory_experiment(mwpm, shots);
        if (world_rank == 0) {
            std::cout << "Expected LER = " << mxp_res.logical_error_rate << "\n";
        }
    }
    sim.load_decoder(mwpm);

    sim.params.verbose = parser.option_set("v") && (world_rank == 0);

    std::string trace_folder;
    if (parser.get_string("trace-folder", trace_folder)) {
        sim.params.save_syndromes_to_file = true;
        sim.params.syndrome_output_folder = trace_folder;
    }

    sim.run(shots);
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

    if (mwpm)   delete mwpm;
}

