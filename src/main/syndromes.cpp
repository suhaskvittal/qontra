/*
 *  author: Suhas Vittal
 *  date:   11 Octoboer 2023
 * */

#include "experiments.h"
#include "instruction.h"
#include "parsing/cmd.h"

#include <deque>
#include <fstream>
#include <iostream>
#include <vector>

#include <mpi.h>

using namespace qontra;
using namespace experiments;

stim::Circuit
get_circuit(const schedule_t& sch, fp_t p) {
    const uint n = get_number_of_qubits(sch);

    tables::ErrorAndTiming et;
    et = et * (1000 * p);
    ErrorTable errors;
    TimeTable timing;
    tables::populate(n, errors, timing, et);
    stim::Circuit circ = schedule_to_stim(sch, errors, timing);
    return circ;
}

int main(int argc, char* argv[]) {
    CmdParser pp(argc, argv);

    std::string asm_file;
    std::string output_file;
    uint64_t shots;

    fp_t p = 1e-3;
    uint64_t hw_min = 1;

    if (!pp.get_string("asm", asm_file)) return 1;
    if (!pp.get_string("out", output_file)) return 1;
    if (!pp.get_uint64("shots", shots)) return 1;

    pp.get_uint64("hw-min", hw_min);
    pp.get_float("p", p);

    G_FILTERING_HAMMING_WEIGHT = hw_min;

    // Get stim circuit.
    schedule_t sch = schedule_from_file(asm_file);
    stim::Circuit circuit = get_circuit(sch, p);

    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Write results to file.
    std::filesystem::path output_path(output_file);
    if (world_rank == 0) {
        std::filesystem::path output_folder(output_path.parent_path());
        safe_create_directory(output_folder);
    }
    std::ofstream out(output_path);
    // Setup experiemnts callback for writing syndromes to file.
    std::deque<std::vector<uint>> syndrome_list;
    callback_t cb;
    cb.prologue = [&] (stim::simd_bits_range_ref& row) {
        auto detectors = get_nonzero_detectors(row, circuit.count_detectors());
        bool nonzero_obs = row.popcnt() - detectors.size();
        if (detectors.size() >= G_FILTERING_HAMMING_WEIGHT || nonzero_obs) {
            syndrome_list.push_back(detectors);
        }
    };
    // Execute experiment.
    generate_syndromes(circuit, shots, cb);
    // Now, each MPI process takes turns writing to the file.
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
