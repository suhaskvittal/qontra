#!/bin/sh
 
# author: Suhas Vittal
# date: 16 November 2023

output_file=../data/hex_color_code_threshold_circuit_level.csv
shots=$2

proc=$1

# Run all of the experiments.
cd build

make -j${proc}

for d in 3
do
    echo "d = ${d}"
    for p in 1e-3
    do
        echo "   p = ${p}"
        # Write memory asm.
        mpirun -np $proc ./memory \
            --asm ../data/asm/memory_x_d${d}.asm \
            --out $output_file \
            --p $p \
            --shots $shots \
            -circuit
    done
done

