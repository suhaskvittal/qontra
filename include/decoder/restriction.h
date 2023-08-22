/* author: Suhas Vittal
 *  date:   19 July 2023
 *
 *  An implementation of Kubica and Delfosse's
 *  Restriction Decoder for color codes.
 * */

#ifndef RESTRICTION_DECODER_h
#define RESTRICTION_DECODER_h

#include "decoder/mwpm.h"
#include "defs.h"
#include "graph/graph.h"
#include "graph/algorithms/distance.h"

#include <bitset>
#include <deque>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace qontra {

#define __COLOR graph::decoding::vertex_t::Color

namespace restriction {

inline int c2i(__COLOR x) {
    if (x == __COLOR::red)      return 0;
    if (x == __COLOR::green)    return 1;
    if (x == __COLOR::blue)     return 2;
    if (x == __COLOR::notred)   return 3;
    if (x == __COLOR::notgreen) return 4;
    if (x == __COLOR::notblue)  return 5;
    return -1;
}

inline __COLOR i2c(int x) {
    if (x == 0) return __COLOR::red;
    if (x == 1) return __COLOR::green;
    if (x == 2) return __COLOR::blue;
    if (x == 3) return __COLOR::notred;
    if (x == 4) return __COLOR::notgreen;
    if (x == 5) return __COLOR::notblue;
    return __COLOR::none;
}

inline __COLOR get_remaining_color(__COLOR c1, __COLOR c2) {
    uint8_t available = 0x7;
    available &= ~(c1 == __COLOR::red || c2 == __COLOR::red);
    available &= ~((c1 == __COLOR::green || c2 == __COLOR::green) << 1);
    available &= ~((c1 == __COLOR::blue || c2 == __COLOR::blue) << 2);
    if (available & 0x1)   return __COLOR::red;
    if (available & 0x2)   return __COLOR::green;
    if (available & 0x4)   return __COLOR::blue;
    return __COLOR::none;
}

typedef std::pair<uint, __COLOR>                    cdet_t;
typedef std::tuple<cdet_t, cdet_t, __COLOR>         match_t;
typedef std::pair<std::vector<match_t>, __COLOR>    component_t;

typedef std::pair<cdet_t, cdet_t>                   cdetpair_t;

inline bool
pair_is_in(std::set<cdetpair_t>& s, cdetpair_t x) {
    cdetpair_t ix = std::make_pair(x.second, x.first);
    return s.count(x) || s.count(ix);
}

inline void
insert_pair_into(std::set<cdetpair_t>& s, cdetpair_t x) {
    cdetpair_t ix = std::make_pair(x.second, x.first);
    s.insert(x);
    s.insert(ix);
}

inline void
xor_pair_into(std::set<cdetpair_t>& s, cdetpair_t x) {
    cdetpair_t ix = std::make_pair(x.second, x.first);
    if (s.count(ix) || s.count(x)) {
        s.erase(ix);
        s.erase(x);
    } else {
        s.insert(x);
        s.insert(ix);
    }
}

// Structure Graph:
struct stv_t : graph::base::vertex_t {
    uint    qubit;
    cdet_t  detector;
    bool    is_flag;
};

struct ste_t : graph::base::edge_t {};

typedef graph::Graph<stv_t, ste_t>  StructureGraph;
typedef graph::distance::DistanceMatrix<stv_t, std::vector<stv_t*>>  StructureMatrix;

}   // restriction

#define CDET_TO_ID(x)   (((x).first) | (((uint64_t)(x).second + 1) << 48))

class RestrictionDecoder : public Decoder {
public:
    RestrictionDecoder(const stim::Circuit&, 
                        const uint detectors_per_round,
                        restriction::StructureGraph*,
                        const std::map<uint64_t, uint64_t>& detector_to_base,
                        const std::map<uint64_t, uint64_t>& flag_to_owner);

    ~RestrictionDecoder(void) {
        for (auto dec : rlatt_dec)  delete dec;
    }

    std::string name(void) override { return "RestrictionDecoder"; }

    Decoder::result_t   decode_error(const syndrome_t&) override;
private:
    // The decode_restricted_lattice function decodes the restricted lattice that
    // restricts the passed in color.
    Decoder::result_t   decode_restricted_lattice(const std::vector<uint>&, __COLOR);

    void    update_weights(const std::vector<uint>& flags);
    void    revert_weights(void);

    bool
        is_flag(restriction::cdet_t);
    std::set<restriction::cdetpair_t>        
        get_all_edges_in_component(const restriction::component_t&);
    std::vector<restriction::cdet_t>
        get_detectors_between(restriction::cdet_t, restriction::cdet_t);
    std::set<restriction::cdet_t>            
        get_common_neighbors(restriction::cdet_t, restriction::cdet_t);
    std::set<restriction::cdet_t>            
        get_neighbors(restriction::cdet_t);
    stim::simd_bits             
        get_correction_for_face(restriction::cdet_t, restriction::cdet_t, restriction::cdet_t);

    restriction::StructureGraph*    lattice_structure;
    restriction::StructureMatrix    lattice_matrix;

    std::array<std::set<restriction::cdet_t>, 3>
                                    boundary_adjacent;

    std::array<MWPMDecoder*, 3> rlatt_dec;  // Three MWPM decoders, each of which responsible
                                            // for decoding a different restricted lattice.
                                            //
                                            // Note RestrictionDecoder is a friend of MWPMDecoder,
                                            // so the decoder can access their DecodingGraphs.
    typedef std::map<uint, std::array<uint, 3>>     fwdetmap_t;
    typedef std::map<std::pair<uint, uint>, uint>   bwdetmap_t;
    // As the MWPM decoders operate on
    // different decoding graphs and these
    // graphs do not match the input
    // circuit, we must be able to convert
    // the input syndrome to syndromes for
    // each MWPM decoder.
    fwdetmap_t  to_rlatt; 
    bwdetmap_t  from_rlatt;

    std::map<uint64_t, __COLOR>     color_map;

    uint                            detectors_per_round;
    std::map<uint64_t, uint64_t>    detector_to_base;   // A map to the base detector, which is
                                                        // the lowest round detector (essentially,
                                                        // detector modulo measurement).
    std::map<uint64_t, uint64_t>    flag_to_owner;  // A map of the flag's base detector to its
                                                    // owning check's base detector.
    std::map<graph::decoding::edge_t*, fp_t>    edge_to_old_weight;
    std::vector<uint>                           affected_decoding_graphs;
};

}   // qontra

#endif  // RESTRICTION_DECODER_h
