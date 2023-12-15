#!/bin/sh

proc=$1
shots=$2

echo "RSC NO LEAKAGE---------------------------------------------"

./scripts/probius/rsc_memory_no_leakage.sh $proc $shots

echo "RSC LEAKAGE, NO LRCS---------------------------------------"

./scripts/probius/rsc_memory_leakage_no_lrcs.sh $proc $shots

echo "RSC ALWAYS LRCS, SWAP--------------------------------------"

./scripts/probius/rsc_memory_always_lrcs_swaps.sh $proc $shots

echo "RSC ALWAYS LRCS, DQLR--------------------------------------"

./scripts/probius/rsc_memory_always_lrcs_dqlr.sh $proc $shots

echo "RSC OPTIMAL LRCS, SWAP-------------------------------------"

./scripts/probius/rsc_memory_optimal_lrcs_swaps.sh $proc $shots

echo "RSC OPTIMAL LRCS, DQLR-------------------------------------"

./scripts/probius/rsc_memory_optimal_lrcs_dqlr.sh $proc $shots
