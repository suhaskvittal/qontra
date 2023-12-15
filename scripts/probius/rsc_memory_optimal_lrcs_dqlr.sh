#!/bin/sh

POLICY_NAME="optimal_lrcs_dqlr"
LEAKAGE_OPTIONS="-enl --lrc-policy optimal -dqlr"

(
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
