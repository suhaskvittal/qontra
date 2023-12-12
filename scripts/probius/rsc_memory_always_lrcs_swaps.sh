#!/bin/sh

POLICY_NAME="always_lrcs_swaps"
LEAKAGE_OPTIONS="-enl --lrc-policy always"

PROC=$1
SHOTS=$2

(
    export PROC
    export SHOTS
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
