/*
 *  author: Suhas Vittal
 *  date:   3 July 2024
 * */

#include <iostream>

#include <vtils/mat2.h>
#include <vtils/timer.h>

using namespace vtils;

int main(int argc, char* argv[]) {
    Mat2 m(4,7);
    // Steane's code:
    m.set(0,0,true);
    m.set(0,1,true);
    m.set(0,2,true);
    m.set(0,3,true);

    m.set(3,2,true);
    m.set(3,3,true);
    m.set(3,4,true);
    m.set(3,5,true);

    m.set(2,1,true);
    m.set(2,2,true);
    m.set(2,5,true);
    m.set(2,6,true);
    
    // And a redundant check for testing:
    m.set(1,0,true);
    m.set(1,3,true);
    m.set(1,5,true);
    m.set(1,6,true);

    Timer tt;
    double t;

    // Compute pivots:
    tt.clk_start();
    auto v1 = get_basis_vectors(m.trr());
    t = tt.clk_end() / 1000.0;
    std::cout << "Pivots ------------------- (time = " << t << ")\n";
    for (auto x : v1) {
        std::cout << x;
    }
    // Compute null basis vectors
    tt.clk_start();
    auto v2 = get_null_basis_vectors(m);
    t = tt.clk_end() / 1000.0;
    std::cout << "Null basis --------------- (time = " << t << ")\n";
    for (auto x : v2) {
        std::cout << x;
    }

    Mat2 m2(7,7);
    for (size_t i = 0; i < v1.size(); i++) {
        for (size_t j = 0; j < 7; j++) {
            m2.set(i,j,v1[i](0,j));
        }
    }
    for (size_t i = 0; i < v2.size(); i++) {
        for (size_t j = 0; j < 7; j++) {
            m2.set(i+v1.size(),j,v2[i](0,j));
        }
    }
    auto v3 = get_basis_vectors(m2.trr());
    std::cout << "Operators --------------- (time = " << t << ")\n";
    for (auto x : v3) {
        std::cout << x;
    }
}
