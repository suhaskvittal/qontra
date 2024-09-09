/*
 *  author: Suhas Vittal
 *  date:   26 August 2024
 * */

namespace cgen {

template <int R, int S, size_t GW, size_t GH, size_t CW, size_t CH>
std::vector< std::tuple<Star<R,S>*, typename MCM::S_eqc_t*, coord_t> >
MCM::get_candidates_for_star_check( Star<R,S>* rt, const coord_t& rt_loc, sample_t& samp) {
    std::vector< std::tuple<Star<R,S>*, S_eqc_t*, coord_t> > cand;

    const uint16_t mllsq = static_cast<uint16_t>( SQR(config.max_link_len) );
    constexpr uint16_t cwsq = SQR(CW), chsq = SQR(CH);
    constexpr int CBWdiv2 = static_cast<int>(CBW) / 2;
    constexpr int CBHdiv2 = static_cast<int>(CBH) / 2;

    size_t best_cycle_size = 0;
    for (S_eqc_t* c_p : samp.star_classes) {
        coord_t base_delta = c_sub(c_p->repr_loc, rt_loc);
        uint16_t init_dsq = SQR(gx16(base_delta)) + SQR(gy16(base_delta));
        if (init_dsq > mllsq) continue;  // No slack.
        uint16_t remaining_slack = mllsq - init_dsq;
        // Add the representative to the candidate list.
        if (c_p->repr != rt) cand.emplace_back(c_p->repr, c_p, PZERO);
        // Now try to move to other cells.
        // Distance is toroidal (for now).
        for (int i = -CBWdiv2; i < CBWdiv2; i++) {
            uint8_t ui = static_cast<uint8_t>( i < 0 ? CBW+i : i );
            if (SQR(i) > remaining_slack) continue;
            for (int j = -CBHdiv2; j < CBHdiv2; j++) {
                if (i == 0 && j == 0) continue;
                uint8_t uj = static_cast<uint8_t>( j < 0 ? CBH+j : j );
                size_t k = ui*CBH + uj;
                if (SQR(i*CW) + SQR(j*CH) <= remaining_slack
                    && k < c_p->elements.size()
                    && c_p->elements[k]->get_size() < R)
                {
                    Star<R,S>* x = c_p->elements.at(k);
                    coord_t off = make_coord(ui, uj);
                    // Compute potential cycle size.
                    Star<R,S>* z = nullptr;
                    size_t cyc = speculate_star_tie_and_get_cycle_size(rt, x, &z, samp);
                    if (cyc >= best_cycle_size) {
                        if (cyc > best_cycle_size) {
                            cand.clear();
                            best_cycle_size = cyc;
                        }
                        cand.emplace_back(x, c_p, off);
                        if (z != nullptr) {
                            samp.star_max_cycle_map[z] = 
                                std::max(samp.star_max_cycle_map[z], cyc);
                        }
                    }
                }
            }
        }
    }
    if (best_cycle_size > 0) {
        samp.min_star_cycle = std::min(samp.min_star_cycle, best_cycle_size);
    }
    return cand;
}

template <int R, int S, size_t GW, size_t GH, size_t CW, size_t CH> size_t
MCM::speculate_star_tie_and_get_cycle_size(
        Star<R,S>* x,
        Star<R,S>* y,
        Star<R,S>** star_with_cyc_p,
        const sample_t& s) 
{
    // y is incoming into x.
    size_t min_cycle = std::numeric_limits<size_t>::max();

    const auto& cte_x = s.cycle_table.at(x),
              & cte_y = s.cycle_table.at(y);
    for (const auto& [z, cnt] : cte_x) {
        if (cte_y.count(z)) {
            size_t cyc = cnt + cte_y.at(z) + 1;
            if (cyc < min_cycle) {
                min_cycle = cyc;
                *star_with_cyc_p = z;
            }
        }
    }
    return min_cycle;
}

template <int R, int S, size_t GW, size_t GH, size_t CW, size_t CH> CycleBuffer<R,S>
MCM::search_for_simple_cycles_upto(
        size_t max_cycle_len,
        Star<R,S>* from,
        const sample_t& samp) 
{
    CycleBuffer<R,S> buf;
    buf.prev[from] = { nullptr, nullptr, 0 };

    std::deque< Star<R,S>* > bfsq{from};
    while (bfsq.size()) {
        Star<R,S>* x = bfsq.front();
        bfsq.pop_front();
        const cycle_prev_entry_t<R,S>& pe_x = buf.prev.at(x);

        for (const sptr<Qubit<R,S>>& q : x->get_qubits_ref()) {
            if (q == nullptr) continue;
            for (Star<R,S>* y : q->x_checks) {
                if (y == nullptr || x == y || y == pe_x.s_ptr) continue;
                if (buf.prev.count(y)) {
                    // This is a cycle. Create a buffer entry.
                    const cycle_prev_entry_t<R,S>& pe_y = buf.prev.at(y);
                    // Compute latest common ancestor.
                    uint8_t xlen = pe_x.len,
                            ylen = pe_y.len;
                    Star<R,S>* xcurr = x,
                             * ycurr = y;
                    // We want to get xcurr and ycurr to reach the same
                    // level in the BFS tree.
                    uint8_t size = 1;
                    size += demote_ptrs_in_bfs_tree(buf, &xcurr, xlen, ylen);
                    size += demote_ptrs_in_bfs_tree(buf, &ycurr, ylen, xlen);
                    // We do not need to track xlen and ylen as xlen = ylen = 0
                    // is the root of the tree.
                    Star<R,S>* anc;
                    while (size < max_cycle_len) {
                        if (xcurr == ycurr) {
                            anc = xcurr;
                            break;
                        } else {
                            // Demote both pointers.
                            for (Star<R,S>* z : {xcurr, ycurr}) {
                                const cycle_prev_entry_t<R,S>& e = buf.prev.at(z);
                                z = e.s_ptr;
                            }
                            size += 2;
                        }
                    }
                    if (size < max_cycle_len) {
                        buf.data.push_back({x, y, anc, q, size});
                    } else {
                        std::cout << "cycle too large: " << size+0 << std::endl;
                    }
                } else {
                    buf.prev[y] = { x, q, pe_x.len+1 };
                    bfsq.push_back(y);
                }
            }
        }
    }
    return buf;
}

}   // cgen
