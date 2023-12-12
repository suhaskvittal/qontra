#!/bin/sh

POLICY_NAME="always_lrcs_dqlr"
LEAKAGE_OPTIONS="-enl --lrc-policy always -dqlr"

PROC=$1
SHOTS=$2

(
    export PROC
    export SHOTS
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
