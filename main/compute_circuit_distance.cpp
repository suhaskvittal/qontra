/*
 *  author: Suhas Vittal
 *  date:   1 April 2024
 * */

#include <qontra/ext/stim.h>
#include <qontra/experiments/memory.h>

using namespace qontra;

int main(int argc, char* argv[]) {
    std::string qes_file(argv[1]);
    DetailedStimCircuit circuit = make_default_circuit(qes_file, 1e-3, true);
    auto dem = 
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
                circuit,
                false,
                true,
                false,
                0.0,
                false,
                false);
    auto errors = stim::find_undetectable_logical_error(dem, 3, 3, false);
    std::cout << "Distance = " << errors.count_errors() << std::endl;
    std::cout << errors << "\n";
}
