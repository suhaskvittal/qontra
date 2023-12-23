#!/bin/sh

POLICY_NAME="always_lrcs_swaps"
LEAKAGE_OPTIONS="-enl --lrc-policy always"

(
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
