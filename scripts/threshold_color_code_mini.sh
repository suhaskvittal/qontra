#!/bin/sh
 
# author: Suhas Vittal
# date: 16 November 2023

output_file=../data/hex_color_code_threshold_mini.csv
shots=$2

proc=$1

# Run all of the experiments.

cd build

make -j${proc}

dpr_m=3
dpr_i=3
dpr_f=6
for d in 5
do
#   dpr=$(( 3*dpr_m ))
    dpr=$(( d*dpr_f ))
    echo "d = ${d} (dpr = ${dpr}, m = ${dpr_m}, i = ${dpr_i})"
    for p in 6e-4 7e-4 8e-4
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
    dpr_f=$(( dpr_f+3 ))
done

