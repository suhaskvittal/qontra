/*
 *  author: Suhas Vittal
 *  date:   2 September 2023
 *
 *  A simple compatability wrapper
 *  for PyMatching.
 * */

#ifndef QONTRA_PYMATCHING_h
#define QONTRA_PYMATCHING_h

#include "qontra/decoder.h"

#include <stim.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>

namespace qontra {

pm::Mwpm init_solver_from_circuit(stim::Circuit);

class PyMatching : public Decoder {
public:
    PyMatching(const DetailedStimCircuit&);

    Decoder::result_t decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) override;
private:
    pm::Mwpm solver;
};

}   // qontra

#endif  // QONTRA_PYMATCHING_h
