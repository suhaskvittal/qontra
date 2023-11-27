#!/bin/sh
 
# author: Suhas Vittal
# date: 16 November 2023

cd scripts
mkdir asm

for d in 3 5 7 9
do
    python asm_gen_color_code.py asm/memory_z_d${d}.asm $d $d 
done
