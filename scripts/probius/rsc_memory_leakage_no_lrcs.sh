#!/bin/sh

POLICY_NAME="leakage_no_lrcs"
LEAKAGE_OPTIONS="-enl"

PROC=$1
SHOTS=$2

(
    export PROC
    export SHOTS
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
