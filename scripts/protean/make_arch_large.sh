#!/bin/sh

declare -a CALLS=(
    "hysc/60_8_6 4"
    "hysc/180_20_9 9"
    "hysc/96_18_8 6"
    "hysc/150_32_6 6"
    "hysc/120_34_5 5"
    "hycc/192_20_8 8 -use-jid -color-checks"
    "hycc/96_20_8 6 -use-jid -color-checks"
)

for ((i = 0; i < ${#CALLS[@]}; i++))
do
    c=${CALLS[$i]}
    echo "scripts/protean/make_arch.py ${c}"
    python3 scripts/protean/make_arch.py ${c}
done

