#!/bin/sh

FAMILY=$1
MIN_SIZE=$2
MAX_SIZE=$3

declare -a HYSC_SF=(
    '4_5' '4_6' '5_5' '5_6'
)

declare -a HYCC_SF=(
    '4_6' '4_8' '4_10' '5_8'
)

echo "${FAMILY} ${MIN_SIZE} ${MAX_SIZE}"

if [[ "$FAMILY" = "hysc" || "$FAMILY" = "all" ]]; then
    for sf in ${HYSC_SF[@]} 
    do
        python scripts/protean/evals.py hysc $sf 2 $MIN_SIZE $MAX_SIZE
    done
fi

if [[ "$FAMILY" = "hycc" || "$FAMILY" = "all" ]]; then
    for sf in ${HYCC_SF[@]} 
    do
        python scripts/protean/evals.py hycc $sf 2 $MIN_SIZE $MAX_SIZE
    done
fi

#if [[ "$FAMILY" = "hexcc" || "$FAMILY" = "all" ]]; then
#    python scripts/protean/evals.py hexcc DNE 0 $MIN_SIZE $MAX_SIZE
#fi
