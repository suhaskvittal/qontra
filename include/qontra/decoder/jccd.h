/*
 *  author: Suhas Vittal
 *  date:   28 September 2024
 *
 *  Joint Color Code Decoder.
 * */

#ifndef QONTRA_DECODER_JCCD_h
#define QONTRA_DECODER_JCCD_h

namespace qontra {

/*
 * Joint Color Code Decoder uses a two-level decoder
 * structure:
 *  (1) First, we execute the `ConcatMWPMDecoder` on
 *  a given syndrome. If the correction weight is
 *  less than a given amount (half the code distance),
 *  then the correction is returned.
 *  (2) However, if the correction is too large, then
 *  we run the second decoder (should be either
 *  `RestrictionDecoder` or `MobiusDecoder`) and
 *  return that correction instead.
 * */
template <class DEC2_T> 
class JointColorCodeDecoder : public Decoder {
public:
    JointColorCodeDecoder(const DetailedStimCircuit& circuit, size_t distance)
        :Decoder(circuit),
        d1(circuit),
        d2(circuit),
        distance(distance)
    {}

    inline Decoder::result_t
    decode_error(stim::simd_bits_range_ref<SIMD_WIDTH> syndrome) override {
        auto r = d1.decode_error(syndrome);
        if (r.correction_weight == distance/2 + 1) {
            auto _r = d2.decode_error(syndrome);
            r.exec_time += _r.exec_time;
            r.corr = std::move(_r.corr);
        }
        return r;
    }
private:
    ConcatMWPMDecoder   d1;
    DEC2_T              d2;

    const size_t distance;
};

}   // qontra

#endif  // QONTRA_DECODER_JCCD_h
