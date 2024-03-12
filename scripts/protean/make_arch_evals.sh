#!/bin/sh

declare -a CALLS=(
# SURFACE CODES -- SMALL
    "hysc/4_5/160_18_8_6 6"
)

for ((i = 0; i < ${#CALLS[@]}; i++))
do
    c=${CALLS[$i]}
    echo "scripts/protean/make_arch.py ${c}"
    python3 scripts/protean/make_arch.py ${c}
done

