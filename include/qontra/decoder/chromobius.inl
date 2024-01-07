/*
 *  author: Suhas Vittal
 *  date:   25 December 2023
 * */

namespace qontra {

inline chromobius::Decoder
init_chromobius(stim::Circuit circuit) {
    const chromobius::DecoderConfigOptions options;

    stim::DetectorErrorModel dem = 
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
            circuit,
            false, // decompose_errors
            true,  // fold loops
            false, // allow gauge detectors
            0.0,   // approx disjoint errors threshold
            false, // ignore decomposition failures
            false
        );
    return chromobius::Decoder::from_dem(dem, options);
}

inline Decoder::result_t
Chromobius::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) {
    timer.clk_start();
    uint64_t _corr = 
        backing_decoder.decode_detection_events({syndrome.u8, syndrome.u8 + syndrome.num_u8_padded()});
    stim::simd_bits<SIMD_WIDTH> corr(1);
    *corr.u64 = _corr;
    fp_t t = (fp_t) timer.clk_end();
    return (Decoder::result_t) {
        t,
        corr
    };
}


}   // qontra
