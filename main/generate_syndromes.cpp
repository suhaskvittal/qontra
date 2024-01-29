/*
 *  author: Suhas Vittal
 *  date:   11 Octoboer 2023
 * */

#include <qontra/experiments.h>
#include <qontra/experiments/memory.h>
#include <qontra/experiments/stats.h>
#include <qontra/ext/qes.h>
#include <qontra/ext/stim.h>

#include <vtils/cmd_parse.h>
#include <vtils/filesystem.h>

#include <deque>
#include <fstream>
#include <iostream>
#include <vector>

#include <mpi.h>

using namespace qontra;

int main(int argc, char* argv[]) {
    vtils::CmdParser pp(argc, argv);

    std::string qes_file;
    std::string output_file;
    uint64_t shots;

    fp_t p = 1e-3;
    uint64_t hw_min = 1;

    pp.get("qes", qes_file, true);
    pp.get("out", output_file, true);
    pp.get("shots", shots, true);

    pp.get("hw-min", hw_min);
    pp.get("p", p);

    G_FILTERING_HAMMING_WEIGHT = hw_min;

    // Get stim circuit.
    DetailedStimCircuit circuit = make_circuit(qes_file, p);

    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Write results to file.
    if (world_rank == 0) {
        vtils::safe_create_directory(vtils::get_parent_directory(output_file.c_str()));
    }
    std::ofstream out(output_file);
    // Generate syndromes and write to file:
    std::deque<std::vector<uint64_t>> syndrome_list;
    statistic_t<fp_t> hw_sum(MPI_SUM), hw_max(MPI_MAX);

    generate_syndromes(circuit, shots,
        [&] (shot_payload_t payload) {
            stim::simd_bits<SIMD_WIDTH> syndrome = payload.syndrome,
                                        observable = payload.observables;
            auto detectors = get_nonzero_detectors_(syndrome, circuit.count_detectors());
            bool nonzero_obs = observable.not_zero();
            if (detectors.size() >= G_FILTERING_HAMMING_WEIGHT || nonzero_obs) {
                syndrome_list.push_back(detectors);
            }
            hw_sum += syndrome.popcnt();
            hw_max.scalar_replace_if_better_extrema(syndrome.popcnt());
        });
    // Now, each MPI process takes turns writing to the file.
    hw_sum.reduce();
    hw_max.reduce();
    statistic_t<fp_t> hw_mean = hw_sum.get_mean(shots);

    if (world_rank == 0) {
        out << "Statistics: mean = " << hw_sum.at() << ", max = " << hw_max.at() << "\n";
    }
    for (int r = 0; r < world_size; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == r) {
            for (const auto& dets : syndrome_list) {
                const size_t hw = dets.size();
                out << "hw = " << hw << "\t|";
                for (uint64_t d : dets) out << "\t" << d;
                out << "\n";
            }
        }
    }
    MPI_Finalize();
    return 0;
}
