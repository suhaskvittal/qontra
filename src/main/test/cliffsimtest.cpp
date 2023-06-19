
#include "sim/clifford_sim.h"

int main(int argc, char* argv[]) {
    using namespace qontra;

    uint64_t s = 4096;
    CliffordSimulator sim(3, s);
    sim.H({0});
    sim.CX({0,1});
    sim.CX({1,2});
    sim.H({2});
    sim.M({0,1,2});

    stim::simd_bit_table meas = sim.record_table.transposed();
    for (uint64_t i = 0; i < s; i++) {
        std::cout << "t" << i << " : " 
            << meas[i][0] 
            << meas[i][1] 
            << meas[i][2] << "\n";
    }
}
