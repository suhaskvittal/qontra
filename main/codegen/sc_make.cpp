/*
 *  author: Suhas Vittal
 *  date:   26 August 2024
 * */

#include <codegen/surface_codes/mc.h>

#include <iostream>

using namespace cgen;

#define R 5
#define S 5
#define GL 6
#define CL 2

#define PERC(x,y) static_cast<fp_t>(x)/static_cast<fp_t>(y)*100.0

int main(int argc, char* argv[]) {
    MonteCarloManager<R,S,GL,GL,CL,CL> mgr;
    mgr.make_sample();
    // Print out performance stats:
    const fp_t perc_time_in_init = PERC(mgr.stats.total_time_in_init, mgr.stats.total_time_overall);

    std::cout << "Performance of SC_MAKE(" << R << ", " << S << ") with grid size "
        << GL << " / " << CL << ": " << "\tTotal time taken: " << mgr.stats.total_time_overall << "ns\n";
    std::cout << "-----------------------------------------------------------------------\n";
    std::cout << "\tTotal time in init:\t" << mgr.stats.total_time_in_init << "ns (" 
                << perc_time_in_init << "%)\n";
    return 0;
}
