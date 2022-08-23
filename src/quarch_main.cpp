/*
 *  author: Suhas Vittal
 *  date:   1 August 2022
 * */

#include "quarch.h"
#include "gulliver/experiments.h"

int main() {
    gulliver_mwpm_timing_experiment();
    gulliver_bfu_timing_experiment();
    gulliver_decoder_analysis_experiment();
    return 0;
}
