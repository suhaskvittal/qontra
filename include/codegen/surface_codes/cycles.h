/*
 *  author: Suhas Vittal
 *  date:   26 August 2024
 * */

#ifndef CGEN_SC_CYCLES_h
#define CGEN_SC_CYCLES_h

#include "codegen/surface_codes/utils.h"

#include <unordered_map>
#include <vector>

namespace cgen {

template <int R, int S> using cycle_table_t=std::unordered_map< Star<R,S>*,  std::unordered_map<Star<R,S>*, size_t> >;
template <int R, int S> using cycle_t=std::array< sptr<Qubit<R,S>>, 8 >;

template <int R, int S>
struct cycle_prev_entry_t {
    Star<R,S>* s_ptr;
    sptr<Qubit<R,S>> qubit;
    uint8_t len;
};

template <int R, int S>
struct cycle_buf_entry_t {
    Star<R,S>* src;
    Star<R,S>* dst;
    Star<R,S>* common_ancestor;
    sptr<Qubit<R,S>> back_qubit;
    uint8_t size;
};

template <int R, int S>
class CycleBuffer {
public:
    cycle_t<R,S> pop_next_buf();

    std::vector<cycle_buf_entry_t<R,S>> data;
    std::unordered_map< Star<R,S>*, cycle_prev_entry_t<R,S> > prev;
private:
    template <bool FWD, class ITER>
    void walk(ITER, Star<R,S>* from, Star<R,S>* anc);
};

// Demotes x_p until xlen == ylen. Returns number of demotions.
template <int R, int S>
uint8_t demote_ptrs_in_bfs_tree(const CycleBuffer<R,S>&, Star<R,S>** x_p, uint8_t xlen, uint8_t ylen);

}   // cgen

#include "cycles.inl"

#endif // CGEN_SC_CYCLES_h
