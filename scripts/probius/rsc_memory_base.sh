#!/bin/sh

if [ "${TRACE}" == "yes" ]; then
    ./scripts/probius/rsc_memory_tracer_base.sh
fi

if [ "${DECODE}" == "yes" ]; then
    ./scripts/probius/rsc_memory_decode_base.sh
fi
