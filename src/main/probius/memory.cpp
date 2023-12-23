/*
 *  author: Suhas Vittal
 *  date:   11 December 2023
 * */

#ifdef QONTRA_PYMATCHING_ENABLED

#include <decoder/pymatching.h>

#else

#include <decoder/mwpm.h>

#endif

#include <experiments.h>
#include <parsing/cmd.h>

#include <fstream>
#include <iostream>

#include <mpi.h>
#include <stdio.h>

using namespace qontra;

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    CmdParser pp(argc, argv);

    std::string error_model_file;
    std::string trace_folder;
    std::string output_file;

    if (!pp.get_string("error-model", error_model_file)) return 1;
    if (!pp.get_string("traces", trace_folder)) return 1;
    if (!pp.get_string("out", output_file)) return 1;

    // Retrieve the stim circuit (error model).
    FILE* error_model_fin = fopen(error_model_file.c_str(), "r");
    DetailedStimCircuit error_model(stim::Circuit::from_file(error_model_fin));
    fclose(error_model_fin);

    // Create the folder containing the output file if it does not exist.
    std::filesystem::path output_file_path(output_file);
    if (world_rank == 0) {
        std::filesystem::path output_folder(output_file_path.parent_path());
        safe_create_directory(output_folder);
    }
    bool write_header = !std::filesystem::exists(output_file_path);

    MPI_Barrier(MPI_COMM_WORLD);
    std::ofstream fout(output_file, std::ios::app);
    if (write_header && world_rank == 0) {
        fout << "Trace,"
                << "Logical Error Rate\n";
    }
    // Run memory experiment.
#ifdef QONTRA_PYMATCHING_ENABLED
    PyMatching dec(error_model);
#else
    MWPMDecoder dec(error_model);
#endif
    experiments::G_FILTER_OUT_SYNDROMES = false;

    experiments::memory_params_t params;
    params.trace_folder = trace_folder;

    auto result = memory_experiment(&dec, params);

    if (world_rank == 0) {
        fout << std::filesystem::path(trace_folder).filename() << ","
            << result.logical_error_rate << "\n";
    }
    MPI_Finalize();
    return 0;
}

