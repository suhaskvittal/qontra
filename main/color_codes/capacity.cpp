/*
 *  author: Suhas Vittal
 *  date:   12 May 2024
 * */

#include "gen.h"

#include <qontra/graph/tanner_graph.h>
#include <qontra/graph/io.h>
#include <qontra/experiments/memory.h>
#include <qontra/ext/stim.h>

#include <qontra/decoder/chromobius.h>
#include <qontra/decoder/restriction.h>
#include <qontra/decoder/neural.h>

#include <vtils/cmd_parse.h>
#include <vtils/filesystem.h>

using namespace qontra;
using namespace graph;

using namespace vtils;

int main(int argc, char* argv[]) {
    configure_optimal_batch_size();
    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::string tanner_graph_file(argv[1]);
    std::string output_file(argv[2]);
    
    CmdParser pp(argc, argv, 2);
    std::string HELP =
        "usage: ./cc_capacity <tanner-graph-file> <output-file>\n"
        "optional:\n"
        "\t--e <errors until stopping, default=25>\n"
        "\t--pmin <error-rate, default=1e-3>\n"
        "\t--pmax <error-rate, default=1e-3>\n"
        "\t--step-size <default=1>\n"
        "\t-x or -z (memory x or memory z, default is \'-z\')\n";
    pp.help = HELP;

    uint64_t errors_until_stop = 25;
    fp_t pmin = 1e-3,
         pmax = 1e-3;
    uint64_t step_size = 1;

    pp.get("e", errors_until_stop);
    pp.get("pmin", pmin);
    pp.get("pmax", pmax);
    pp.get("step-size", step_size);

    bool mx = pp.option_set("x");

    // Get tanner graph.
    TannerGraph tgr = 
        create_graph_from_file<TannerGraph>(tanner_graph_file, io::update_tanner_graph);
    std::map<sptr<tanner::vertex_t>, int> color_map;
    if (tgr.compute_check_color_map(color_map) > 2) {
        std::cerr << "Found >3-coloring." << std::endl;
        exit(1);
    }
    fp_t p = pmin;
    DetailedStimCircuit base_circuit = make_capacity(&tgr, p, mx, color_map);
    RestrictionDecoder dec(base_circuit);

    memory_config_t config;
    config.errors_until_stop = errors_until_stop;

    while (p <= 1.1*pmax) {
        DetailedStimCircuit circuit = make_capacity(&tgr, p, mx, color_map);
        dec.set_circuit(circuit);
        memory_result_t res = run_memory_with_generated_syndromes(&dec, config);
        // Write result to file.
        if (world_rank == 0) {
            std::cout << "Writing p = " << p << std::endl;
            bool write_header = false;
            if (!file_exists(get_parent_directory(output_file.c_str()))) {
                safe_create_directory(get_parent_directory(output_file.c_str()));
                write_header = true;
            }
            std::ofstream fout(output_file, std::ios::app);
            if (write_header) {
                fout << "tanner graph,physical error rate,block error rate,"
                    "norm. block error rate" << std::endl;
            }
            fout << get_basename(tanner_graph_file) << ","
                << p << ","
                << res.logical_error_rate << ","
                << res.logical_error_rate / static_cast<fp_t>(base_circuit.count_observables());
            fout << std::endl;
        }
        fp_t gran = floor(log10(p));
        p += step_size*pow(10, gran);
    }
    MPI_Finalize();
    return 0;
}
