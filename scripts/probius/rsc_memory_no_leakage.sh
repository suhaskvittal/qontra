#!/bin/sh

POLICY_NAME="no_leakage"
LEAKAGE_OPTIONS=""

PROC=$1
SHOTS=$2

(
    export PROC
    export SHOTS
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
