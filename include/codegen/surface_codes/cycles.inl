/*
 *  author: Suhas Vittal
 *  date:   26 August 2024
 * */

namespace cgen {

template <int R, int S> cycle_t<R,S>
CycleBuffer<R,S>::pop_next_buf() {
    cycle_buf_entry_t<R,S> e = std::move(data.back());
    data.pop_back();

    cycle_t<R,S> cyc;

    walk<true>(cyc.begin(), e.src, e.common_ancestor);
    walk<false>(cyc.begin()+e.size-1, e.dst, e.common_ancestor);
    cyc[e.size-1] = e.back_qubit;
    return cyc;
}

template <int R, int S>
template <bool FWD, class ITER> inline void
CycleBuffer<R,S>::walk(ITER it, Star<R,S>* curr, Star<R,S>* anc) {
    constexpr int delta = FWD ? 1 : -1;

    while (curr != anc) {
        const cycle_prev_entry_t<R,S>& e = prev.at(curr);
        *it = e.qubit;
        curr = e.s_ptr;
        std::advance(it, delta);
    }
    const cycle_prev_entry_t<R,S>& e = prev.at(anc);
    *it = e.qubit;
}

template <int R, int S> inline uint8_t
demote_ptrs_in_bfs_tree(
        const CycleBuffer<R,S>& buf,
        Star<R,S>** x_p,
        uint8_t xlen,
        uint8_t ylen)
{
    const uint8_t old_xlen = xlen;
    while (xlen > ylen) {
        const cycle_prev_entry_t<R,S>& e = buf.prev.at(*x_p);
        *x_p = e.s_ptr;
        xlen--;
    }
    return old_xlen-xlen;
}

} // cgen
