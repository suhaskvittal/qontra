#!/bin/sh

declare -a CALLS=(
# SURFACE CODES -- SMALL
#   "hysc/4_5/160_18_8_6 6"
# COLOR CODES -- SMALL
    "hycc/4_6/24_8_4_4 4 -use-jid -color-checks"
    "hycc/4_6/96_20_6_6 6 -use-jid -color-checks"
#   "hycc/4_6/216_40_8_8 8 -use-jid -color-checks"
)

for ((i = 0; i < ${#CALLS[@]}; i++))
do
    c=${CALLS[$i]}
    echo "scripts/protean/make_arch.py ${c}"
    python3 scripts/protean/make_arch.py ${c}
done

