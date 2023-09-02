/*
 *  author: Suhas Vittal
 *  date:   1 September 2023
 * */

#ifndef ASTREA_DECODER_h
#define ASTREA_DECODER_h

#include "decoder/decoder.h"
#include "graph/decoding_graph.h"

namespace qontra {
namespace mldh {

template <class D_t, int K=10>
class AstreaDecoder : public Decoder {
public:
    AstreaDecoder(const stim::Circuit& circuit)
        :Decoder(circuit, graph::DecodingGraph::Mode::LOW_MEMORY),
        global_decoder(circuit)
    {}

    Decoder::result_t decode_error(const syndrome_t& syndrome) override {
        auto res = global_decoder.decode_error(syndrome);
        
        const uint hw = get_nonzero_detectors(syndrome).size();
        
        fp_t t;
        if (hw <= K)        t = 4;
        else if (hw <= K+2) t = 4;  // Can be parallelized.
        else if (hw <= K+4) t = (K+3)*4;
        else if (hw <= K+6) t = ((K+5)*(K+3))*4;
        else                t = res.exec_time;

        return (Decoder::result_t) {
            t,
            res.corr,
            res.is_error,
            res.error_assignments
        };
    }

    std::string name(void) { return "Astrea_K" + std::to_string(K); }
private:
    D_t global_decoder;
};

}   // mldh
}   // qontra

#endif  // ASTREA_DECODER_h
