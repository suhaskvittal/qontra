#!/bin/sh

POLICY_NAME="always_lrcs_dqlr"
LEAKAGE_OPTIONS="-enl --lrc-policy always -dqlr"

(
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
