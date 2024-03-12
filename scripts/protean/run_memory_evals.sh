#!/bin/sh

nproc=$1

declare -a CALLS=(
# SURFACE CODES -- SMALL
    "hysc/4_5/160_18_8_6 mwpm"
)

for ((i = 0; i < ${#CALLS[@]}; i++))
do
    c=${CALLS[$i]}
    memcmd="scripts/protean/run_memory.py ${c} 5e-4 5e-3 10000 ${nproc} -fix-error"
    echo $memcmd
    python3 $memcmd
done

