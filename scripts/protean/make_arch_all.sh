#!/bin/sh

declare -a CALLS=(
    "hysc/30_8_3 3"
    "hysc/36_8_4 4"
    "hysc/48_6_4 4"
    "hysc/60_18_4 4"
    "hycc/24_8_4 4 -use-jid -color-checks"
    "hycc/48_8_8 4 -use-jid -color-checks"
)

for ((i = 0; i < ${#CALLS[@]}; i++))
do
    c=${CALLS[$i]}
    echo "scripts/protean/make_arch.py ${c}"
    python3 scripts/protean/make_arch.py ${c}
done
