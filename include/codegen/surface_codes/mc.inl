/*
 *  author: Suhas Vittal
 *  date:   26 August 2024
 * */

namespace cgen {

#define MCM MonteCarloManager<R,S,GW,GH,CW,CH>

#define GET_RAND_IT(arr) arr.begin() + (rng()%arr.size())
#define SQR(x) (x)*(x)

template <int R, int S, size_t W, size_t H> void
StarClass<R,S,W,H>::connect_with_delta(
        StarClass<R,S,W,H>* other,
        const coord_t& delta,
        uint8_t qi,
        uint8_t qj) 
{
    uint16_t off = dim<static_cast<uint16_t>(H)>(delta);
    for (size_t i = 0; i < TOTAL_COUNT; i++) {
        Star<R,S>* a = elements.at(i),
                 * b = other->elements.at( (i+off)%TOTAL_COUNT );
        a->tie(b, qi, qj);
    }
}

template <int R, int S, size_t GW, size_t GH, size_t CW, size_t CH> MCM::sample_t
MCM::make_sample() {
    global_timer.clk_start();

    sample_t s = init_sample();

    stats.total_time_overall += global_timer.clk_end();
    return s;
}

template <int R, int S, size_t GW, size_t GH, size_t CW, size_t CH> MCM::sample_t
MCM::init_sample() {
    local_timer.clk_start();

    sample_t samp;
    // Initialize the sample with the appropriate parameters.
    for (size_t i = 0; i < CW; i++) {
        for (size_t j = 0; j < CH; j++) {
            if (fpdst(rng) > config.star_init_prob) {
                continue;
            }
            Star<R,S>* repr = new Star<R,S>;
            S_eqc_t* c_p = new S_eqc_t(repr, make_coord(i,j));
            for (Star<R,S>* x : c_p->elements) {
                samp.star_checks.push_back(x);
                samp.s_eqc_map[x] = c_p;
            }
            samp.star_classes.push_back(std::move(c_p));
        }
    }
    // Add links between star classes.
    for (S_eqc_t* c_p : samp.star_classes) {
        Star<R,S>* ch1 = c_p->repr;
        auto cand = get_candidates_for_star_check(ch1, c_p->repr_loc, samp);
        for (size_t j = ch1->get_size(); j < R; j++) {
            if (cand.empty()) break;
            auto r_it = GET_RAND_IT(cand);
            // Connect the two checks.
            auto& [ch2, ch2_c_p, delta] = *r_it;
            // Compute delta between ch1 and ch2.
            c_p->connect_with_delta(ch2_c_p, delta, j, ch2->get_size());
            cand.erase(r_it);
        }
    }

    stats.total_time_in_init += local_timer.clk_end();
    return samp;
}

template <int R, int S, size_t GW, size_t GH, size_t CW, size_t CH>
std::vector< std::tuple<Star<R,S>*, typename MCM::S_eqc_t*, coord_t> >
MCM::get_candidates_for_star_check( Star<R,S>* rt, const coord_t& rt_loc, const sample_t& samp) {
    std::vector< std::tuple<Star<R,S>*, S_eqc_t*, coord_t> > cand;

    const uint16_t mllsq = static_cast<uint16_t>( SQR(config.max_link_len) );
    constexpr uint16_t cwsq = SQR(CW), chsq = SQR(CH);
    constexpr int CBWdiv2 = static_cast<int>(CBW) / 2;
    constexpr int CBHdiv2 = static_cast<int>(CBH) / 2;
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
                if (SQR(i) + SQR(j) <= remaining_slack
                    && k < c_p->elements.size())
                {
                    coord_t off = make_coord(ui, uj);
                    cand.emplace_back(c_p->elements.at(k), c_p, off);
                }
            }
        }
    }
    return cand;
}

}   // cgen
