/*
 *  author: Suhas Vittal
 *  date:   12 February 2024
 * */

#include <qontra/experiments/memory.h>
#include <qontra/graph/decoding_graph.h>
#include <qontra/ext/stim.h>
#include <qontra/ext/qes.h>

#include <iomanip>
#include <iostream>

#include <stdlib.h>

using namespace qontra;

// This file is a potpourri file. Modify as you'd like.

int main(int argc, char* argv[]) {
    std::string qes_file(argv[1]);
    fp_t p = atof(argv[2]);

    DetailedStimCircuit error_model = make_circuit(qes_file, p, true);

    uptr<graph::DecodingGraph> gr = std::make_unique<graph::DecodingGraph>(error_model, 3);
    /*
    // Print out error polynomial:
    poly_t errp = gr.get_error_polynomial();

    std::cout << " ====== Error Polynomial ====== " << std::endl;;
    for (size_t i = 0; i < errp.size() && i < 50; i++) {
        std::cout << std::setw(6) << i << " | "
            << std::setw(8)
            << std::setprecision(6)
            << std::scientific
            << errp.at(i)
            << std::endl;
    }
    std::cout << " ------------------------------ " << std::endl;;
    std::cout << std::setw(6) << "exp" << " | "
        << std::setw(8)
        << std::setprecision(6)
        << gr.get_expected_errors()
        << std::endl;
    */
    
    return 0;
}
