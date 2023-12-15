#!/bin/sh

POLICY_NAME="no_leakage"
LEAKAGE_OPTIONS=""

(
    export POLICY_NAME
    export LEAKAGE_OPTIONS

    ./scripts/probius/rsc_memory_base.sh
)
