#!/bin/sh
 
# author: Suhas Vittal
# date: 16 November 2023

output_file=../data/hex_color_code_threshold_phenomenological.csv
shots=$2

proc=$1

# Run all of the experiments.
cd build

make -j${proc}

dpr_m=1
dpr_i=2
for d in 3
do
    dpr=$(( 3*dpr_m ))
    echo "dpr = ${dpr}, m = ${dpr_m}, i = ${dpr_i}"
    for p in 0.0 8e-3 1e-2 2e-2
    do
        # Write memory asm.
        mpirun -np $proc ./memory \
            --asm ../data/asm/memory_z_d${d}.asm \
            --out $output_file \
            --p $p \
            --shots $shots \
            --dpr $dpr \
            -pheno
    done
    dpr_m=$(( dpr_m+dpr_i ))
    dpr_i=$(( dpr_i+1 ))
done

