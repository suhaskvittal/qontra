#!/bin/sh
 
# author: Suhas Vittal
# date: 16 November 2023

output_file=../data/hex_color_code_threshold_mini.csv
shots=$2

proc=$1

# Run all of the experiments.

cd build

make -j8

dpr_m=3
dpr_i=3
for d in 5 7
do
    dpr=$(( 3*dpr_m ))
    echo "d = ${d} (dpr = ${dpr}, m = ${dpr_m}, i = ${dpr_i})"
    for p in 6e-3 7e-3 8e-3
    do
        echo "  p = ${p}"
        # Write memory asm.
        mpirun -np $proc ./memory \
            --asm ../scripts/asm/memory_z_d${d}.asm \
            --out $output_file \
            --p $p \
            --shots $shots \
            --dpr $dpr
    done
    dpr_m=$(( dpr_m+dpr_i ))
    dpr_i=$(( dpr_i+1 ))
done

