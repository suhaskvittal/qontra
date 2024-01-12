/*
 *  author: Suhas Vittal
 *  date:   11 January 2024
 * */

#include <qontra/ext/qes.h>

qontra::DetailedStimCircuit
make_circuit(std::string qes_file, fp_t p) {
    using namespace qontra;

    qes::Program<> program = qes::from_file(qes_file);
    const size_t n = get_number_of_qubits(program);

    tables::ErrorAndTiming et;
    et.e_g1q *= 0.1;
    et.e_idle *= 0.1;
    et = et * (1000*p);

    ErrorTable errors;
    TimeTable timing;
    tables::populate(n, errors, timing, et);
    return DetailedStimCircuit::from_qes(program, errors, timing);
}
