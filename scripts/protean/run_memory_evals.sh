#!/bin/sh

nproc=$1

declare -a CALLS=(
# SURFACE CODES -- SMALL
#   "hysc/4_5/60_8_4_6 mwpm"
    "hysc/5_5/30_8_3_3 mwpm"
    "hysc/5_5/80_18_5_5 mwpm"
    "hysc/5_6/60_18_3_4 mwpm"
# COLOR CODES -- SMALL
    "hycc/4_6/24_8_4_4 restriction"
    "hycc/4_6/96_20_7_7 restriction"
    "hycc/4_8/32_12_4_4 restriction"
# SURFACE CODES -- MEDIUM
    "hysc/4_5/160_18_6_8 restriction"
    "hysc/5_6/120_34_5_6 restriction"
# COLOR CODES -- MEDIUM
    "hycc/4_6/192_36_8_8 restriction"
    "hycc/4_8/144_40_8_8 restriction"
)

for ((i = 0; i < ${#CALLS[@]}; i++))
do
    c=${CALLS[$i]}
    memcmd="scripts/protean/run_memory.py ${c} 5e-4 5e-3 10000 ${nproc} -fix-error"
    echo $memcmd
    python3 $memcmd
done

