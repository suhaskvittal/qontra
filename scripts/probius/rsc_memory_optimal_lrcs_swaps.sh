#!/bin/sh

POLICY_NAME="optimal_lrcs_swaps"
LEAKAGE_OPTIONS="-enl --lrc-policy optimal"

(
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
