/*
 *  author: Suhas Vittal
 *  date:   15 July 2023
 *
 *  Testing for a perceptron round predictor.
 * */

#include "graph/io.h"
#include "graph/lattice_graph.h"
#include "experiments.h"
#include "parsing/cmd.h"
#include "sim/components/syndrome_predictor.h"
#include "sim/components/syndrome_predictors/multi_level.h"
#include "sim/components/syndrome_predictors/perceptron.h"

#include <fstream>
#include <iostream>
#include <random>

#include <stim.h>

#include <math.h>

using namespace qontra;
using namespace experiments;
using namespace sim;

int main(int argc, char* argv[]) {
    CmdParser parser(argc, argv);

    G_USE_MPI = false;

    uint distance;
    uint rounds;
    if (!parser.get_uint32("distance", distance)
        || !parser.get_uint32("rounds", rounds))    return 1;
    stim::CircuitGenParameters params(rounds, distance, "rotated_memory_z");
    stim::Circuit circuit = stim::generate_surface_code_circuit(params).circuit;

    std::string trace_folder;
    std::string lattice_file;
    if (!parser.get_string("trace-folder", trace_folder))   return 1;
    if (!parser.get_string("lattice-file", lattice_file))   return 1;

    int detectors_per_round = (distance*distance-1) / 2;

    std::mt19937_64 rng(0);
    std::vector<stim::simd_bits> testing_data;

    int d;
    if (!parser.get_int32("depth", d)) d = 2;
    
    // Create predictors.
    std::ifstream fin(lattice_file);
    graph::io::callback_t<graph::LatticeGraph> cb = 
        [] (graph::LatticeGraph& graph, std::string line)
        {
            graph::io::update_lattice_graph(graph, line);
        };
    graph::LatticeGraph lattice = graph::create_graph_from_file(fin, cb);

    SimpleSyndromePredictor simple(d, detectors_per_round);
    PerceptronPredictor pshare(d,
                                    detectors_per_round,
                                    1,
                                    &lattice,
                                    PerceptronPredictor::Mode::share);
    PerceptronPredictor pgeom(d,
                                    detectors_per_round,
                                    1,
                                    &lattice,
                                    PerceptronPredictor::Mode::geom);
    PerceptronPredictor plocal(d,
                                    detectors_per_round,
                                    1,
                                    &lattice,
                                    PerceptronPredictor::Mode::local);
    PerceptronPredictor pglobal(d,
                                    detectors_per_round,
                                    1,
                                    &lattice,
                                    PerceptronPredictor::Mode::global);

    PerceptronPredictor ml_per(plocal);
    ml_per.predict_only_nontrivial_events = true;

    uint64_t max_data = 50000;

    stim::simd_bit_table data(max_data+100, (d+1)*detectors_per_round);
    stim::simd_bit_table test(max_data+100, (d+1)*detectors_per_round);
    uint64_t s1 = 0;
    uint64_t s2 = 0;

    callback_t callbacks;
    callbacks.prologue = [&] (stim::simd_bits_range_ref row) 
    {
        uint obs_bits = 0;
        for (uint i = 0; i < circuit.count_observables(); i++) {
            obs_bits += row[i+circuit.count_detectors()];
        }
        const uint hw = row.popcnt() - obs_bits;
        if (hw == 0)    return;
        if (rng() % 100 < 20 && s2 < max_data) {
            // This is testing data.
            test[s2++] |= row;
            return;
        }
        if (s1 >= max_data)  return;
        for (uint r = d-1; r < rounds-1; r++) {
            // Create the training data set.
            uint min = (r - (d-1))*detectors_per_round;
            uint max = (r+2)*detectors_per_round;
            bool all_zeros = true;
            for (uint i = 0; i < max-min; i++) {
                all_zeros &= !row[i+min];
            } 
            if (all_zeros)  continue;
            for (uint i = 0; i < max-min; i++) {
                data[s1][i] = row[i+min];
            } 
            s1++;
        }
    };
    read_syndrome_trace(trace_folder, circuit, callbacks);
    data = data.transposed();

    pshare.train(data);
    ml_per.train(data);

    MultiLevelPredictor<2> multi_level;
    multi_level[1] = &simple;
    multi_level[0] = &ml_per;

    /*

    std::cout << "\tgeom...\n";
    pgeom.train(data);

    std::cout << "\tlocal...\n";
    plocal.train(data);

    std::cout << "\tglobal...\n";
    pglobal.train(data);


    std::cout << "Done training.\n";
    */

    std::vector<SyndromePredictor*> predictors;
    predictors.push_back(&simple);
    predictors.push_back(&pshare);
    predictors.push_back(&multi_level);
    /*
    predictors.push_back(&pglobal);
    predictors.push_back(&pgeom);
    predictors.push_back(&plocal);
    */
    pshare.print_weights();

    std::vector<std::string> pnames{
        "simple",
        "perceptron_share",
        "multi_level",
        "perceptron_global",
        "perceptron_geom",
        "perceptron_local",
    };

    auto test_tp = test.transposed();

    std::ofstream   pp_pred_out("predictions.txt");

    const uint off1 = (d-1)*detectors_per_round;
    const uint off2 = d*detectors_per_round;

    // Write to log files before computing stats.
    auto pp_res = multi_level.predict(test_tp);
    auto pp_sig_pred = pp_res.sig_pred.transposed();
    auto pp_val_pred = pp_res.val_pred.transposed();
    for (uint64_t t = 0; t < s2; t++) {
        pp_pred_out << "----------------------------------------------\n";
        for (uint r = 0; r <= d; r++) {
            if (r == d) pp_pred_out << "****";
            else        pp_pred_out << "    ";
            pp_pred_out << "\t\t";
            for (uint i = 0; i < detectors_per_round; i++) {
                if (test[t][r*detectors_per_round+i])   pp_pred_out << "!";
                else                                    pp_pred_out << "_";
            }
            pp_pred_out << "\n";
        }
        pp_pred_out << "perceptron\t";
        for (uint i = 0; i < detectors_per_round; i++) {
            if (!pp_sig_pred[t][i]) {
                pp_pred_out << "?";
            } else if (pp_val_pred[t][i] == test[t][off2+i]) {
                pp_pred_out << "O";
            } else {
                pp_pred_out << "X";
            }
        }
        pp_pred_out << "\n";
    }

    // Compute stats for each predictor.
    for (uint k = 0; k < predictors.size(); k++) {
        auto predictor = predictors[k];
        std::string name = pnames[k];
        auto res = predictor->predict(test_tp);
        
        auto sig_pred = res.sig_pred.transposed();
        auto val_pred = res.val_pred.transposed();

        uint64_t n_predicts = 0;
        uint64_t n_correct_predicts = 0;
        uint64_t n_nontrivial_predicts = 0;
        uint64_t n_correct_nontrivial_predicts = 0;
        uint64_t n_requests = 0;
        for (uint64_t t = 0; t < s2; t++) {
            fp_t n_predicted = 0;
            for (uint i = 0; i < detectors_per_round; i++) {
                n_requests++;
                if (!sig_pred[t][i])    continue;
                n_predicts++;
                n_correct_predicts += (val_pred[t][i] == test[t][off2+i]);
                n_nontrivial_predicts += test[t][off1+i];
                n_correct_nontrivial_predicts +=
                        (val_pred[t][i] == test[t][off2+i])
                        && (test[t][off1+i]);
                n_predicted++;
            }
        }

        std::cout << "[ " << name << " ]\n";
        std::cout << "\tPredicted " << n_predicts << " of " 
                << n_requests << " (" 
                << (((fp_t)n_predicts)/n_requests * 100.0)
                << "%)\n";
        std::cout << "\tCorrectly predicted " << n_correct_predicts
                << " of " << n_predicts << " ("
                << (((fp_t)n_correct_predicts)/n_predicts * 100.0)
                << "%)\n";
        std::cout << "\t\t(nontrivial) " << n_correct_nontrivial_predicts
                << " of " << n_nontrivial_predicts << " ("
                << (((fp_t)n_correct_nontrivial_predicts)/n_nontrivial_predicts
                        * 100.0) << "%)\n";
    }
}
