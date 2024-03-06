#!/bin/sh

nproc=$1

declare -a CALLS=(
    "hysc/60_8_6 mwpm 1e-3 1e-3 100"
#   "hysc/180_20_9 mwpm 1e-3 6e-2 10000"
#   "hysc/96_18_8 mwpm 5e-4 3e-3 10000"
#   "hysc/150_32_6 mwpm 5e-4 3e-3 10000"
#   "hysc/120_34_5 mwpm 5e-4 3e-3 10000"
#   "hycc/192_20_8 restriction 4e-4 1e-3 10000"
#   "hycc/96_20_8 restriction 4e-4 1e-3 10000"
)

for ((i = 0; i < ${#CALLS[@]}; i++))
do
    c=${CALLS[$i]}
    echo "scripts/protean/run_memory.py ${c} ${nproc} -fix-error"
    python3 scripts/protean/run_memory.py ${c} ${nproc} -fix-error
done

