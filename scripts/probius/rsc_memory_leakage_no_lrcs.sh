#!/bin/sh

POLICY_NAME="leakage_no_lrcs"
LEAKAGE_OPTIONS="-enl"

(
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
