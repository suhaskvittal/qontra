#!/bin/sh
#
# author: Suhas Vittal
# date: 16 November 2023

output_file=../data/hex_color_code_threshold.csv
shots=$2

proc=$1

# Write all of the ASM files.
conda deactivate
conda activate venv

cd scripts
mkdir asm

for d in 3 5 7 9
do
    python asm_gen_color_code.py asm/memory_z_d${d}.asm $d $d 
done

conda deactivate

cd ../build

dpr_m=1
dpr_i=2
for d in 3 5 7 9
do
    dpr=$(( 3*dpr_m ))
    echo "dpr = ${dpr}, m = ${dpr_m}, i = ${dpr_i}"
    for p in 5e-4 6e-4 7e-4 8e-4 9e-4 1e-3 2e-3 3e-3 4e-3 5e-3
    do
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

