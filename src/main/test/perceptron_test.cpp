/*
 *  author: Suhas Vittal
 *  date:   15 July 2023
 *
 *  Testing for a perceptron round predictor.
 * */

#include "experiments.h"
#include "parsing/cmd.h"

#include <fstream>
#include <iostream>
#include <random>

#include <stim.h>

#include <math.h>

using namespace qontra;
using namespace experiments;

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
    if (!parser.get_string("trace-folder", trace_folder))   return 1;

    int detectors_per_round = (distance*distance-1) / 2;

    std::mt19937_64 rng(0);
    std::vector<stim::simd_bits> testing_data;

    int d;
    if (!parser.get_int32("depth", d)) d = 2;
    
    //
    //  Two perceptrons so far:
    //      (1) All-zeros prediction.
    //      (2) Same-bits prediction.

    std::vector<int32_t>    weights_allzeros(detectors_per_round*d+1, 0);
    std::vector<int32_t>    weights_issame(detectors_per_round*d+1, 0);

    std::ofstream fout("syndrome_log.txt");

    callback_t callbacks;
    callbacks.prologue = [&] (stim::simd_bits_range_ref row) 
    {
        uint obs_bits = 0;
        for (uint i = 0; i < circuit.count_observables(); i++) {
            obs_bits += row[i+circuit.count_detectors()];
        }
        const uint hw = row.popcnt() - obs_bits;
        if (hw == 0)    return;
        fout << "--------------------------------------------------------\n";
        if (rng() % 100 < 20) {
            // This is testing data.
            stim::simd_bits row_cpy(row);
            testing_data.push_back(row_cpy);
            return;
        }
        for (int r = 0; r < rounds; r++) {
            fout << "R" << r << "\t";
            for (uint i = 0; i < detectors_per_round; i++) {
                uint k = r*detectors_per_round + i;
                if (row[k]) fout << "!";
                else        fout << "_";
            }
            fout << "\n";
        }
        for (int r = d-1; r < rounds-1; r++) {
            // We collect the true labels for the training data.
            bool is_all_zeros = true;
            bool is_all_same = true;
            for (uint i = 0; i < detectors_per_round; i++) {
                is_all_zeros &= !row[(r+1)*detectors_per_round + i];
                is_all_same &= row[r*detectors_per_round+i]
                                    == row[(r+1)*detectors_per_round+i];
            }
            // Update weights accordingly.
            int y_allzeros = is_all_zeros ? 1 : -1;
            int y_issame = is_all_same ? 1 : -1;
            weights_allzeros[0] += y_allzeros;
            weights_issame[0] += y_issame;
            int wn = 1;
            for (int i = -(d-1)*detectors_per_round; i < detectors_per_round; i++) {
                int rv = row[r*detectors_per_round + i] ? 1 : -1;
                weights_allzeros[wn] += y_allzeros * rv;
                weights_issame[wn] += y_issame * rv;
                wn++;
            }
        }
    };
    read_syndrome_trace(trace_folder, circuit, callbacks);

    uint correct_pred_allzeros = 0;
    uint total_pred_allzeros = 0;
    uint correct_pred_naive = 0;
    uint correct_nonzero_pred_allzeros = 0;
    uint total_nonzero_pred_allzeros = 0;

    uint correct_pred_issame = 0;
    uint total_pred_issame = 0;
    uint correct_nonzero_pred_issame = 0;
    uint total_nonzero_pred_issame = 0;

    for (auto row : testing_data) {
        for (int r = d-1; r < rounds-1; r++) {
            bool is_all_zeros = true;
            bool is_all_same = true;
            bool curr_round_is_all_zeros = true;
            for (uint i = 0; i < detectors_per_round; i++) {
                is_all_zeros &= !row[(r+1)*detectors_per_round + i];
                is_all_same &= row[r*detectors_per_round+i]
                                    == row[(r+1)*detectors_per_round + i];
                curr_round_is_all_zeros &= !row[r*detectors_per_round+i];
            }
            int true_y_allzeros = is_all_zeros ? 1 : -1;
            int true_y_issame = is_all_same ? 1 : -1;
            int pred_y_allzeros = weights_allzeros[0];
            int pred_y_issame = weights_issame[0];
            int wn = 1;
            bool input_is_all_zeros = true;
            for (int i = -(d-1)*detectors_per_round; i < detectors_per_round; i++) {
                int rv = row[r*detectors_per_round + i] ? 1 : -1;
                pred_y_allzeros += weights_allzeros[wn] * rv;
                pred_y_issame += weights_issame[wn] * rv;

                input_is_all_zeros &= (rv < 0);

                wn++;
            }
            pred_y_allzeros = pred_y_allzeros > 0 ? 1 : -1;
            pred_y_issame = pred_y_issame > 0 ? 1 : -1;
            
            correct_pred_allzeros += true_y_allzeros == pred_y_allzeros;
            correct_pred_naive += is_all_zeros && curr_round_is_all_zeros;
            total_pred_allzeros++;
            correct_nonzero_pred_allzeros += 
                !curr_round_is_all_zeros && (true_y_allzeros == pred_y_allzeros);
            total_nonzero_pred_allzeros += !curr_round_is_all_zeros;

            correct_pred_issame += true_y_issame == pred_y_issame;
            total_pred_issame++;
            correct_nonzero_pred_issame += 
                !curr_round_is_all_zeros && (true_y_issame == pred_y_issame);
            total_nonzero_pred_issame += !curr_round_is_all_zeros;
        }
    }
    std::cout << "Perceptron Statistics (All Zeros):\n";
    std::cout << "\t" << correct_pred_allzeros << " of " << total_pred_allzeros
                << " (" 
                << (((fp_t)correct_pred_allzeros)/total_pred_allzeros * 100.0) 
                << "%)\n";
    std::cout << "\t\tnaive approach predicts " << correct_pred_naive
                << " of " << total_pred_allzeros << " (" 
                << (((fp_t)correct_pred_naive)/total_pred_allzeros * 100.0)
                << "%)\n";
    std::cout << "\t(only nonzero) " << correct_nonzero_pred_allzeros << " of " 
                << total_nonzero_pred_allzeros << " (" 
                << (((fp_t)correct_nonzero_pred_allzeros)/total_nonzero_pred_allzeros 
                        * 100.0) 
                << "%)\n";

    std::cout << "Perceptron Statistics (All Same):\n";
    std::cout << "\t" << correct_pred_issame << " of " << total_pred_issame
                << " (" 
                << (((fp_t)correct_pred_issame)/total_pred_issame * 100.0) 
                << "%)\n";
    std::cout << "\t(only nonzero) " << correct_nonzero_pred_issame << " of " 
                << total_nonzero_pred_issame << " (" 
                << (((fp_t)correct_nonzero_pred_issame)/total_nonzero_pred_issame 
                        * 100.0) 
                << "%)\n";
}
