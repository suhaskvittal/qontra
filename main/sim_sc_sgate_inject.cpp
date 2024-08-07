/*
 *  author: Suhas Vittal
 *  date:   18 July 2024
 * */

#include <qontra/sim/base/clifford_sim.h>
#include <qontra/experiments.h>

#include <vtils/utility.h>

#include <mpi.h>
#include <stdlib.h>

using namespace qontra;

inline std::vector<fp_t>
uni_noise(uint64_t n, fp_t p) {
    return std::vector<fp_t>(n, p);
}

bool matches(
    stim::simd_bit_table<SIMD_WIDTH>& tab1,
    stim::simd_bit_table<SIMD_WIDTH>& tab2,
    uint64_t t)
{
    auto t1tr = tab1.transposed();
    auto t2tr = tab2.transposed();
    return t1tr[t] == t2tr[t];
}

void
run_batch(
        uint64_t d,
        fp_t p,
        uint64_t trials,
        const std::vector<uint64_t>& all_qubits,
        const std::vector<uint64_t>& all_data,
        const std::vector<uint64_t>& all_parity,
        const std::vector<uint64_t>& all_x_parity,
        const std::map<uint64_t, std::vector<int64_t>>& cx_order_map,
        const std::set<uint64_t>& xs,
        const std::set<uint64_t>& presat,
        std::vector<CliffordSimulator>& sims,
        uint64_t& postsel_trials,
        uint64_t& good_trials,
        uint64_t& final_trials)
{
    auto& s1 = sims[0],
        & s2 = sims[1];
    for (auto& s : sims) {
        s.shots = trials;
        s.reset_sim();
    }

    // Simulate the state injection.
    const uint64_t qinj = d*(d/2) + d/2;
    std::vector<uint64_t> init_h_ops;
    for (uint64_t r = 0; r < d; r++) {
        for (uint64_t c = 0; c < d; c++) {
            uint64_t q = d*r + c;
            if (r == d/2 && c == d/2
                || r < d/2 && c >= d/2
                || r > d/2 && c <= d/2) 
            {
                init_h_ops.push_back(q);
            }
        }
    }
    for (auto& s : sims) s.H(init_h_ops);
    // Inject depolarizing error on all qubits (state-preparation).
//  s2.error_channel<&StateSimulator::eDP1>( all_data, uni_noise(all_data.size(), 0.1*p) );
    for (auto& s : sims) s.S({qinj});

    // Perform two rounds of syndrome extraction.
    uint64_t mctr = 0;
    std::map<uint64_t, uint64_t> prev_mctr_map, mctr_map;
    for (int r = 0; r < 2; r++) {
        prev_mctr_map = std::move(mctr_map);
        if (r == 1) {
            // Inject decoherence + dephasing error.
//          s2.error_channel<&StateSimulator::eDP1>( all_data, uni_noise(all_data.size(), p) );
        }
        // Reset and init parity qubits.
        for (auto& s : sims) s.R(all_parity);
//      s2.error_channel<&StateSimulator::eX>( all_parity, uni_noise(all_parity.size(), 0.1*p) );
        for (auto& s : sims) s.H(all_x_parity);
//      s2.error_channel<&StateSimulator::eDP1>( all_x_parity, uni_noise(all_x_parity.size(), 0.1*p) );
        // Perform CNOTs.
        for (int k = 0; k < 4; k++) {
            std::vector<uint64_t> cxo;
            for (const auto& [q1,arr] : cx_order_map) {
                int64_t q2 = arr.at(k);
                if (q2 == -1) continue;
                if (xs.count(q1)) {
                    vtils::push_back_all(cxo, {q1, static_cast<uint64_t>(q2)} );
                } else {
                    vtils::push_back_all(cxo, {static_cast<uint64_t>(q2), q1} );
                }
            }
            for (auto& s : sims) s.CX(cxo);
//          s2.error_channel<&StateSimulator::eDP2>( cxo, uni_noise( cxo.size()/2, p ) );
        }
        for (auto& s : sims) s.H(all_x_parity);
        s2.error_channel<&StateSimulator::eDP1>( all_x_parity, uni_noise(all_x_parity.size(), 0.1*p) );

        s1.M(all_parity,
                uni_noise(all_parity.size(), 0), 
                uni_noise(all_parity.size(), 0), 
                mctr);
        s2.M(all_parity,
                uni_noise(all_parity.size(), 0), 
                uni_noise(all_parity.size(), 0), 
                mctr);
        for (uint64_t x : all_parity) {
            mctr_map[x] = mctr++;
        }
    }
    // Reset all parity qubits.
    for (auto& s : sims) s.R(all_parity);
    // First, check if trial is post-selected.
    for (uint64_t tt = 0; tt < trials; tt++) {
        bool is_ps = true;
        // Check post-selection.
        for (uint64_t x : presat) {
            uint64_t m = prev_mctr_map[x];
            if ( s2.record_table[m][tt] ) is_ps = false;
        }
        for (uint64_t x : all_parity) {
            uint64_t m1 = prev_mctr_map[x],
                    m2 = mctr_map[x];
            if ( s2.record_table[m1][tt] ^ s2.record_table[m2][tt] ) is_ps = false;
        }
        bool is_ok = matches(s1.x_table, s2.x_table, tt)
                        && matches(s1.z_table, s2.z_table, tt)
                        && matches(s1.r_table, s2.r_table, tt);
        if (is_ps) postsel_trials++;
        if (is_ok) good_trials++;
        if (is_ps && is_ok) final_trials++;
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(NULL, NULL);
    int world_rank = 0,
        world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    configure_optimal_batch_size();

    uint64_t d = atoi(argv[1]);
    uint64_t t = atoi(argv[2]);
    fp_t p = atof(argv[3]);

    uint64_t lt = t / world_size + (world_rank==0)*(t % world_size);

    uint64_t n = 2*d*d-1;
    std::vector<uint64_t> all_qubits;
    std::vector<uint64_t> all_data;
    std::vector<uint64_t> all_parity;
    std::map<uint64_t, std::vector<int64_t>> cx_order_map;
    std::set<uint64_t> xs;
    // Setup data structures:
    for (uint64_t q = 0; q < d*d; q++) {
        all_qubits.push_back(q);
        all_data.push_back(q);
    }
    uint64_t qs = d*d;
    for (uint64_t r = 0; r < d-1; r++) {
        for (uint64_t c = 0; c < d-1; c++) {
            int64_t q1 = d*r + c,
                    q2 = d*r + c+1,
                    q3 = d*(r+1) + c,
                    q4 = d*(r+1) + c+1;
            if ( (r+c) % 2 == 0 ) {
                cx_order_map[qs] = { q1, q2, q3, q4 };
                xs.insert(qs);
            } else {
                cx_order_map[qs] = { q1, q3, q2, q4 };
            }
            all_parity.push_back(qs++);
        }
    }
    for (uint64_t r = 0; r < d-1; r += 2) {
        int64_t q1 = d*r, 
                q2 = d*(r+1),
                q3 = d*(r+1) + d-1,
                q4 = d*(r+2) + d-1;
        cx_order_map[qs] = { -1, -1, q1, q2 };
        all_parity.push_back(qs++);

        cx_order_map[qs] = { q3, q4, -1, -1 };
        all_parity.push_back(qs++);
    }
    for (uint64_t c = 0; c < d-1; c += 2) {
        int64_t q1 = d*(d-1) + c,
                q2 = d*(d-1) + c+1,
                q3 = c+1,
                q4 = c+2;
        cx_order_map[qs] = { q1, q2, -1, -1 };
        all_parity.push_back(qs);
        xs.insert(qs++);

        cx_order_map[qs] = { -1, -1, q3, q4 };
        all_parity.push_back(qs);
        xs.insert(qs++);
    }
    std::vector<uint64_t> all_x_parity( xs.begin(), xs.end() );
    // Determine presat checks.
    std::set<uint64_t> presat;
    for (const auto& [x, arr] : cx_order_map) {
        bool all_in_same_quadrant = true;
        for (int64_t q : arr) {
            if ( q < 0 ) continue;
            int64_t r = q/d,
                    c = q%d;
            all_in_same_quadrant &= (r != d/2 || c != d/2);
            if (xs.count(x)) {
                all_in_same_quadrant &= (r < d/2 && c >= d/2) || (r > d/2 && c <= d/2);
            } else {
                all_in_same_quadrant &= (r <= d/2 && c < d/2) || (r >= d/2 && c > d/2);
            }
        }
        if (all_in_same_quadrant) presat.insert(x);
    }

    // Do two simulations: one noise free, and one with noise.
    std::vector<CliffordSimulator> sims(2, CliffordSimulator(n,lt));
    uint64_t _t = lt;
    uint64_t _postsel_trials = 0,
             _good_trials = 0,
             _final_trials = 0;
    while (_t) {
        uint64_t trials_this_batch = _t < G_SHOTS_PER_BATCH ? _t : G_SHOTS_PER_BATCH;
        run_batch(d,
                p,
                trials_this_batch,
                all_qubits,
                all_data,
                all_parity,
                all_x_parity,
                cx_order_map,
                xs,
                presat,
                sims,
                _postsel_trials,
                _good_trials,
                _final_trials);
        _t -= trials_this_batch;
    }
    uint64_t postsel_trials,
             good_trials,
             final_trials;
    MPI_Reduce(&_postsel_trials, &postsel_trials, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&_good_trials, &good_trials, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&_final_trials, &final_trials, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        // Report final stats.
        fp_t e = 1.0 - static_cast<fp_t>(final_trials) / static_cast<fp_t>(postsel_trials);
        std::cout << "Post-selected: " << postsel_trials << " of " << t << std::endl;
        std::cout << "Good: " << good_trials << " of " << t << std::endl;
        std::cout << "Good and Post-selected: " << final_trials << " of " << t << std::endl;
        std::cout << "Error rate: " << e << std::endl;
    }
    MPI_Finalize();
    return 0;
}
