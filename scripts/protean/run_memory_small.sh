#!/bin/sh

nproc=$1

declare -a CALLS=(
    "hysc/30_8_3 mwpm 5e-4 3e-3 10000"
    "hysc/36_8_4 mwpm 5e-4 3e-3 10000" 
    "hysc/48_6_4 mwpm 1e-3 1e-2 10000"
    "hysc/60_18_4 mwpm 5e-4 3e-3 10000"
    "hycc/24_8_4 restriction 1e-4 1e-3 100000"
    "hycc/48_8_8 restriction 1e-4 1e-3 100000"
)

for ((i = 0; i < ${#CALLS[@]}; i++))
do
    c=${CALLS[$i]}
    echo "scripts/protean/run_memory.py ${c} ${nproc} -fix-error"
    python3 scripts/protean/run_memory.py ${c} ${nproc} -fix-error
done
