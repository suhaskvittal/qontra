/*
 *  author: Suhas Vittal
 *  date: 2 September 2023
 * */

#include "mldh/latency_sim.h"
#include "parsing/cmd.h"

using namespace qontra;
using namespace mldh;

int main(int argc, char* argv[]) {
    CmdParser pp(argc, argv);

    MPI_Init(NULL, NULL);

    std::string data_folder;
    std::string output_file;
    uint64_t round_latency = 1000;

    if (!pp.get_string("in", data_folder))  return 1;
    if (!pp.get_string("out", output_file)) return 1;
    pp.get_uint64("t", round_latency);

    latency_sim_params_t params;
    params.round_latency = round_latency;

    auto stats = simulate_on_latency_data(data_folder, params);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if (world_rank == 0) {
        std::filesystem::path output_path(output_file);
        std::filesystem::path output_folder(output_path.parent_path());
        safe_create_directory(output_folder);
        bool write_header = !std::filesystem::exists(output_path);
        std::ofstream out(output_path, std::ios::app);
        if (write_header) {
            out << "Data Folder,"
                << "Round Latency,"
                << "Normalized Latency\n";
        }
        out << std::filesystem::path(data_folder).filename() << ","
            << round_latency << ","
            << stats.normalized_latency << "\n";
    }
    MPI_Finalize();
    return 0;
}
