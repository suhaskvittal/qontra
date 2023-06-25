/*
 *  author: Suhas Vittal
 *  date:   23 June 2023
 * */

#include "decoder/mwpm.h"
#include "sim/control_sim.h"

#include <stim.h>

#include <iostream>
#include <vector>

using namespace qontra;

int main(int argc, char* argv[]) {
    schedule_t  prog;
    
    stim::CircuitGenParameters circ_params(3, 3, "rotated_memory_z");
    auto circ = stim::generate_surface_code_circuit(circ_params).circuit;

    const uint n = from_stim_circuit(circ, prog);
    prog.push_back({"decode", {}, {}});
    prog.push_back({"savem", {}, {}});
    prog.push_back({"done", {}, {}});
    /*
    const uint n = 2;
    prog.push_back({"h", {0}, {}});
    prog.push_back({"cx", {0, 1}, {}});
    prog.push_back({"mrc", {0, 1}, {}});
    prog.push_back({"obs", {0, 1}, {}});
    prog.push_back({"obs", {1, 2}, {}});
    prog.push_back({"savem", {}, {}});
    prog.push_back({"done", {}, {}});
    */

    std::cout << schedule_to_text(prog) << "\n";

    experiments::G_USE_MPI = false;
    experiments::G_SHOTS_PER_BATCH = 128;

    ControlSimulator sim(n, prog);
    tables::ErrorAndTiming params;
    tables::populate(n, sim.params.errors, sim.params.timing, params);
    sim.build_canonical_circuit();

    stim::Circuit error_model = sim.get_canonical_circuit();
    std::cout << error_model.str() << "\n";
    decoder::MWPMDecoder mwpm(error_model);

    sim.load_decoder(&mwpm);

    sim.params.verbose = 0;

    sim.run(1L << 16);
    std::cout << "Probability histogram:\n";
    for (auto pair : sim.prob_histograms) {
        std::cout << "\t";
        for (uint64_t x : pair.first)   std::cout << x << " ";
        std::cout << "\t" << pair.second << "\n";
    }
    std::cout << "Time taken to simulate: " << sim.sim_time * 1e-9 << "\n";
}

