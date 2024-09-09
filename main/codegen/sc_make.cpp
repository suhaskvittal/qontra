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
    for (size_t i = 0; i < 1000; i++) mgr.make_sample();
    auto s = mgr.make_sample();
    FILE* fout = fopen("sc_make.out", "w");
    mgr.dump_sample_to_file(fout, s);
    fclose(fout);

    const auto& stats = mgr.stats;
    
    std::cout << "Min star cycle: " << s.min_star_cycle << std::endl;
    std::cout << "Empty samples: " << stats.empty_samples << std::endl;
    std::cout << "No cycles found: " << stats.no_cycles_found << std::endl;

    // Print out performance stats:
    const fp_t perc_time_in_init = PERC(stats.total_time_in_init, stats.total_time_overall);
    const fp_t perc_time_in_cycle_comp = PERC(stats.total_time_in_cycle_comp, stats.total_time_overall);

    std::cout << "Performance of SC_MAKE(" << R << ", " << S << ") with grid size "
        << GL << " / " << CL << ": " << "\tTotal time taken: " << stats.total_time_overall*1e-3 << "us\n";
    std::cout << "-----------------------------------------------------------------------\n";
    std::cout << "\tTotal time in init:\t" << stats.total_time_in_init*1e-3 << "us (" 
                << perc_time_in_init << "%)\n";
    std::cout << "\tTotal time in cycle computation:\t" << stats.total_time_in_cycle_comp*1e-3 << "us (" 
                << perc_time_in_cycle_comp << "%)\n";
    return 0;
}
