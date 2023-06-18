
#include "sim/clifford_sim.h"

int main(int argc, char* argv[]) {
    using namespace qontra;

    uint64_t s = 32;
    CliffordSimulator sim(2, s);
    sim.H({0});
    sim.CX({0,1});
    sim.M({0,1});

    stim::simd_bit_table meas = sim.record_table.transposed();
    for (uint64_t i = 0; i < s; i++) {
        std::cout << "t" << i << " : " << meas[i][0] << meas[i][1] << "\n";
    }
}
