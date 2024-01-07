/*
 *  author: Suhas Vittal
 *  date:   11 Octoboer 2023
 * */

#include <qontra/experiments.h>
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
using namespace experiments;

DetailedStimCircuit
get_circuit(const qes::Program<>& program, fp_t p) {
    const size_t n = get_number_of_qubits(program);

    tables::ErrorAndTiming et;
    et = et * (1000 * p);
    ErrorTable errors;
    TimeTable timing;
    tables::populate(n, errors, timing, et);
    return DetailedStimCircuit::from_qes(program, errors, timing);
}

int main(int argc, char* argv[]) {
    vtils::CmdParser pp(argc, argv);

    std::string qes_file;
    std::string output_file;
    uint64_t shots;

    fp_t p = 1e-3;
    uint64_t hw_min = 1;

    if (!pp.get_string("qes", qes_file)) return 1;
    if (!pp.get_string("out", output_file)) return 1;
    if (!pp.get_uint64("shots", shots)) return 1;

    pp.get_uint64("hw-min", hw_min);
    pp.get_float("p", p);

    G_FILTERING_HAMMING_WEIGHT = hw_min;

    // Get stim circuit.
    qes::Program<> program = qes::from_file(qes_file);
    DetailedStimCircuit circuit = get_circuit(program, p);

    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Write results to file.
    if (world_rank == 0) {
        vtils::safe_create_directory(vtils::get_parent_directory(output_file.c_str()));
    }
    std::ofstream out(output_file);
    // Setup experiemnts callback for writing syndromes to file.
    std::deque<std::vector<uint>> syndrome_list;
    callback_t cb;
    fp_t hw_sum, __hw_sum = 0.0;
    fp_t hw_max, __hw_max = 0.0;
    cb.prologue = [&] (shot_payload_t payload) {
        stim::simd_bits<SIMD_WIDTH> syndrome = payload.syndrome,
                                    observable = payload.observables;

        auto detectors = get_nonzero_detectors_(syndrome, circuit.count_detectors());
        bool nonzero_obs = observable.not_zero();
        if (detectors.size() >= G_FILTERING_HAMMING_WEIGHT || nonzero_obs) {
            syndrome_list.push_back(detectors);
        }
        __hw_sum += syndrome.popcnt();
        if (detectors.size() > __hw_max) __hw_max = detectors.size();
    };
    // Execute experiment.
    generate_syndromes(circuit, shots, cb);
    // Now, each MPI process takes turns writing to the file.
    if (world_size == 1) {
        hw_sum = __hw_sum;
        hw_max = __hw_max;
    }
    MPI_Reduce(&__hw_sum, &hw_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&__hw_max, &hw_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    hw_sum /= shots;
    if (world_rank == 0) {
        out << "Statistics: mean = " << hw_sum << ", max = " << hw_max << "\n";
    }
    for (int r = 0; r < world_size; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (world_rank == r) {
            for (const auto& dets : syndrome_list) {
                const uint hw = dets.size();
                out << "hw = " << hw << "\t|";
                for (auto d : dets) out << "\t" << d;
                out << "\n";
            }
        }
    }
    MPI_Finalize();
    return 0;
}
