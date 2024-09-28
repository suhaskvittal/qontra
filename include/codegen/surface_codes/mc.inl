/*
 *  author: Suhas Vittal
 *  date:   26 August 2024
 * */

#include <algorithm>
#include <deque>
#include <limits>
#include <iostream>

namespace cgen {

#define MCM MonteCarloManager<R,S,GW,GH,CW,CH>

#define GET_RAND_IT(arr) arr.begin() + (rng()%arr.size())
#define SQR(x) (x)*(x)

template <int R, int S, size_t W, size_t H> void
StarClass<R,S,W,H>::connect_with_delta(
        StarClass<R,S,W,H>* other,
        const coord_t& delta,
        cycle_table_t<R,S>& cycle_table) 
{
#ifdef SC_DEBUG
    std::cout << "Connecting C" << print_coord(repr_loc) << " to C" << print_coord(other->repr_loc)
        << " with delta " << print_coord(delta) << "\n";
#endif
    uint16_t off = dim<static_cast<uint16_t>(H)>(delta);
    std::array<size_t, TOTAL_COUNT> visited;
    visited.fill(std::numeric_limits<size_t>::max());
    for (size_t i = 0; i < TOTAL_COUNT; i++) {
        size_t j = (i+off) % TOTAL_COUNT;
        if (this == other && visited[j] == i) break;  // Do not double count.
        Star<R,S>* a = elements.at(i),
                 * b = other->elements.at(j);
        a->tie(b);
        visited[i] = j;
        
        const auto& cte_a = cycle_table.at(a);
        for (const auto& [x, cnt] : cte_a) {
            if (!cycle_table[b].count(x)) {
                cycle_table[b][x] = cnt+1;
            }
        }
    }
}

template <int R, int S, size_t GW, size_t GH, size_t CW, size_t CH> MCM::sample_t
MCM::make_sample() {
    global_timer.clk_start();

    sample_t s = init_sample();
    if (s.star_checks.empty()) {
        stats.empty_samples++;
        goto ms_end;
    }
    add_face_checks(s);

ms_end:
    stats.total_time_overall += global_timer.clk_end();
    return s;
}

template <int R, int S, size_t GW, size_t GH, size_t CW, size_t CH> void
MCM::dump_sample_to_file(FILE* fout, const sample_t& samp) {
    std::unordered_map<sptr<Qubit<R,S>>, size_t> qubit_enum_map;
    std::unordered_map<Star<R,S>*, size_t> star_enum_map;

    size_t n = 0;
    for (Star<R,S>* x : samp.star_checks) {
        star_enum_map[x] = n++;
        // Enumerate the qubits.
        for (size_t i = 0; i < x->get_size(); i++) {
            sptr<Qubit<R,S>> q = x->at(i);
            if (!qubit_enum_map.count(q)) {
                qubit_enum_map[q] = n++;
            }
        }
    }
    // First, write S_eqc_t data.
    for (size_t i = 0; i < samp.star_classes.size(); i++) {
        S_eqc_t* c_p = samp.star_classes[i];
        fprintf(fout, "c%d,%d,%d", i, gx16(c_p->repr_loc), gy16(c_p->repr_loc));
        for (Star<R,S>* x : c_p->elements) {
            fprintf(fout, ",%d", star_enum_map[x]);
        }
        fprintf(fout, "\n");
    }
    // Add connections.
    for (const auto& [q,qi] : qubit_enum_map) {
        fprintf(fout, "d%d", qi);
        for (size_t i = 0; i < q->x_ptr; i++) {
            fprintf(fout, ",%d", star_enum_map[q->x_checks[i]]);
        }
        fprintf(fout, "\n");
    }
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
                // Initialize other structures.
                samp.cycle_table[x][x] = 0;
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
            c_p->connect_with_delta(ch2_c_p, delta, samp.cycle_table);
            cand.erase(r_it);
        }
    }

    stats.total_time_in_init += local_timer.clk_end();
    return samp;
}

template <int R, int S, size_t GW, size_t GH, size_t CW, size_t CH> void
MCM::add_face_checks(sample_t& samp) {
    local_timer.clk_start();

    CycleBuffer<R,S> buf = search_for_simple_cycles_upto(8, samp.star_checks[0], samp);
#ifdef SC_DEBUG
    std::cout << "[ MCM::add_face_checks ] number of cycle buffer entries: " << buf.data.size() << std::endl;
#endif
    if (buf.data.empty()) {
        stats.no_cycles_found++;
        return;
    }
    // Peek the next cycle buffer entry.
    std::cout << "Cycle size: " << buf.data.back().size+0 << "\n";
    cycle_t<R,S> cyc = buf.pop_next_buf();

    stats.total_time_in_cycle_comp += local_timer.clk_end();
}

}   // cgen
