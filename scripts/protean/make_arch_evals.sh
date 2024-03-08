#!/bin/sh

declare -a CALLS=(
# SURFACE CODES -- SMALL
    "hysc/4_5/60_8_4_6 6"
    "hysc/5_5/30_8_3_3 3"
    "hysc/5_5/80_18_5_5 5"
    "hysc/5_6/60_18_3_4 3"
# COLOR CODES -- SMALL
    "hycc/4_6/24_8_4_4 4 -use-jid -color-checks"
    "hycc/4_6/96_20_7_7 7 -use-jid -color-checks"
    "hycc/4_8/32_12_4_4 4 -use-jid -color-checks"
# SURFACE CODES -- MEDIUM
    "hysc/4_5/160_18_6_8 8"
    "hysc/5_6/120_34_5_6 6"
# COLOR CODES -- MEDIUM
    "hycc/4_6/192_36_8_8 8 -use-jid -color-checks"
    "hycc/4_8/144_40_8_8 8 -use-jid -color-checks"
)

for ((i = 0; i < ${#CALLS[@]}; i++))
do
    c=${CALLS[$i]}
    echo "scripts/protean/make_arch.py ${c}"
    python3 scripts/protean/make_arch.py ${c}
done

