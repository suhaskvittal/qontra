/*
 *  author: Suhas Vittal
 *  date:   27 November 2023
 * */

#ifndef QONTRA_DECODER_UTILS_h
#define QONTRA_DECODER_UTILS_h

namespace qontra {

// Below, we define some helper functions for common
// subroutines during decoding.

fp_t    get_edge_weight_between(DecodingGraph&, uint, uint, stim::simd_bits_range_ref);

}   // qontra

#endif  // QONTRA_DECODER_FLAG_h
