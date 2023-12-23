#!/bin/sh

POLICY_NAME="eraser_dqlr"
LEAKAGE_OPTIONS="-enl --lrc-policy eraser -dqlr"

(
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
