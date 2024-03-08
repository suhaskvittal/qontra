/*
 *  author: Suhas Vittal
 *  date:   21 January 2024
 * */

#include <qontra/decoder/mwpm.h>
#include <qontra/decoder/restriction.h>
#include <qontra/ext/stim.h>

#include <qontra/protean/experiments.h>
#include <qontra/experiments.h>
#include <qontra/experiments/memory.h>

#include <vtils/cmd_parse.h>
#include <vtils/filesystem.h>

#ifdef PROTEAN_PERF
#include <vtils/timer.h>
#endif

#include <fstream>
#include <iostream>

#include <mpi.h>

using namespace qontra;
using namespace protean;
using namespace vtils;

int main(int argc, char* argv[]) {
    configure_optimal_batch_size();
    // Initialize MPI
    MPI_Init(NULL, NULL);
    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    CmdParser pp(argc, argv, 1);

    std::string protean_folder(argv[1]);
    std::string qes_file = pp.option_set("mx") ? 
                                protean_folder + "/memory/xrm1.qes" : protean_folder + "/memory/zrm1.qes";
    std::string coupling_file = protean_folder + "/coupling_graph.txt";
    std::string output_file = pp.option_set("mx") ?
                                protean_folder + "/output/basic_memory_x.csv"
                                : protean_folder + "/output/basic_memory_z.csv";

    std::string decoder_name;

    uint64_t    shots = 1'000'000;
    fp_t        pmin = 5e-4,
                pmax = 3e-3;
    uint64_t    step_size = 1;

    pp.get("decoder", decoder_name, true);

    pp.get("s", shots);
    pp.get("pmin", pmin);
    pp.get("pmax", pmax);
    pp.get("step-size", step_size);

    G_BASE_SEED = 1000;

#ifdef PROTEAN_PERF
    Timer timer;
    fp_t t;

    timer.clk_start();
#endif

    // Initialize error and timing tables.
    qes::Program<> program = qes::from_file(qes_file);

#ifdef PROTEAN_PERF
    t = timer.clk_end();
    if (world_rank == 0) {
        std::cout << "[ pr_base_memory ] took " << t*1e-9 << "s to read the QES file." << std::endl;
    }
    timer.clk_start();
#endif

    ErrorTable errors;
    TimeTable timing;
    make_error_and_timing_from_coupling_graph(coupling_file, errors, timing);

    DetailedStimCircuit base_circuit;
    if (pp.option_set("fix-error")) {
        base_circuit = make_circuit(program, pmax);
    } else {

        ErrorTable errors_base = errors * 1e-3;
        TimeTable timing_base = timing * 1e-3;
        base_circuit = DetailedStimCircuit::from_qes(program, errors_base, timing_base);
    }

    uptr<Decoder> dec = nullptr;
    if (decoder_name == "mwpm") {
        dec = std::make_unique<MWPMDecoder>(base_circuit);
    } else if (decoder_name == "restriction") {
        dec = std::make_unique<RestrictionDecoder>(base_circuit);
    } else {
        std::cerr << "Unsupported decoder type: " << decoder_name << std::endl;
    }

#ifdef PROTEAN_PERF
    t = timer.clk_end();
    if (world_rank == 0) {
        std::cout << "[ pr_base_memory ] took " << t*1e-9 << "s to initialize the decoder." << std::endl;
    }
#endif

    fp_t p = pmin;
    while (p <= 1.1*pmax) {
        DetailedStimCircuit circuit;
        // Load model from file and run memory experiment.
        if (pp.option_set("fix-error")) {
            circuit = make_circuit(program, p);
        } else {
            ErrorTable _errors = errors * p;
            TimeTable _timing = timing * p;
            circuit = DetailedStimCircuit::from_qes(program, _errors, _timing);
        }
        dec->set_circuit(circuit);

        memory_config_t config;
        config.shots = shots;
#ifdef PROTEAN_PERF
        timer.clk_start();
#endif
        auto res = memory_experiment(dec.get(), config);
#ifdef PROTEAN_PERF
        t = timer.clk_end();
        std::cout << "[ pr_base_memory ] took " << t*1e-9 << "s to decode at p = " << p << std::endl;
#endif

        // Write result to file.
        if (world_rank == 0) {
            std::cout << "[ pr_base_memory ] Writing p = " << p << std::endl;
            bool write_header = !file_exists(output_file);
            if (!file_exists(get_parent_directory(output_file.c_str()))) {
                safe_create_directory(get_parent_directory(output_file.c_str()));
            }
            std::ofstream fout(output_file, std::ios::app);
            if (write_header) {
                fout << "physical error rate,shots,logical error rate,word error rate" << std::endl;
            }
            fout << p << ","
                << shots << ","
                << res.logical_error_rate << ","
                << res.logical_error_rate / static_cast<fp_t>(circuit.count_observables());
            fout << std::endl;
        }
        fp_t gran = floor(log10(p));
        p += step_size*pow(10, gran);
    }
    MPI_Finalize();
    return 0;
}

