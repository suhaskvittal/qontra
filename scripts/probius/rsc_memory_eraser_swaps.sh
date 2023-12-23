#!/bin/sh

POLICY_NAME="eraser_swaps"
LEAKAGE_OPTIONS="-enl --lrc-policy eraser"

(
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
