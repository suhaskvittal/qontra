/*
 *  author: Suhas Vittal
 *  date:   27 November 2023
 * */

#ifndef QONTRA_EXT_STIM_h
#define QONTRA_EXT_STIM_h

#include "qontra/tables.h"

#include <stim.h>
#include <qes.h>

#include <map>
#include <set>

namespace qontra {

// We provide an extension of Stim's Circuit class. The DetailedStimCircuit
// class can be used in any function accepting stim::Circuits, but has
// extra metadata that can be used by decoders and such.
//
// Stim's circuits cannot hold much information beyond error rates
// and operations, so we contain extra information in here.
struct DetailedStimCircuit : public stim::Circuit {
    DetailedStimCircuit();
    DetailedStimCircuit(const stim::Circuit&);
    DetailedStimCircuit(const DetailedStimCircuit&);
    DetailedStimCircuit(DetailedStimCircuit&&);

    DetailedStimCircuit& operator=(const stim::Circuit&);
    DetailedStimCircuit& operator=(const DetailedStimCircuit&);

    // Set fix_timing_error_as_depolarizing_error to something >= 0 to use depolarizing errors with
    // error rates equal to that value. Otherwise, timing errors are implemented via the Pauli
    // twirling approximation for amplitude and phase damping.
    static DetailedStimCircuit from_qes(const qes::Program<>&, const ErrorTable&, const TimeTable&, 
                                        fp_t fix_timing_error_as_depolarizing_error=-1.0);

    void    apply_errors(std::string, std::vector<uint64_t>, std::vector<fp_t>, bool is_2q);

    int                             number_of_colors_in_circuit;
    std::map<uint64_t, uint64_t>    detector_base_map;
    std::map<uint64_t, int>         detector_color_map;
    std::set<uint64_t>              flag_detectors;
    std::map<uint64_t, uint64_t>    flag_owner_map;
};

// Extension of stim::simd_bits_range_ref for_each_word function.
template <size_t W, size_t N, class FUNC>
void for_each_word(std::array<stim::simd_bits_range_ref<W>, N>, FUNC);

// Template and specialization of ANDNOT for general use.
bool                                    andnot(bool x, bool y);
template <size_t W> stim::bitword<W>    andnot(stim::bitword<W>, stim::bitword<W>);

template <size_t W> void    left_shift(stim::simd_bit_table<W>&, int64_t by);
template <size_t W> void    right_shift(stim::simd_bit_table<W>&, int64_t by);

}   // qontra

namespace stim {

template <size_t W>
bitword<W> operator~(bitword<W>);
    
}   // stim

#include "stim.inl"

#endif  // QONTRA_EXT_STIM_h
