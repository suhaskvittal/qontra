#!/bin/sh

source scripts/probius_setup.sh $1 $2

echo "RSC NO LEAKAGE---------------------------------------------"

./scripts/probius/rsc_memory_no_leakage.sh

echo "RSC LEAKAGE, NO LRCS---------------------------------------"

./scripts/probius/rsc_memory_leakage_no_lrcs.sh

echo "RSC ALWAYS LRCS, SWAP--------------------------------------"

./scripts/probius/rsc_memory_always_lrcs_swaps.sh

echo "RSC ALWAYS LRCS, DQLR--------------------------------------"

./scripts/probius/rsc_memory_always_lrcs_dqlr.sh

echo "RSC OPTIMAL LRCS, SWAP-------------------------------------"

./scripts/probius/rsc_memory_optimal_lrcs_swaps.sh

echo "RSC OPTIMAL LRCS, DQLR-------------------------------------"

./scripts/probius/rsc_memory_optimal_lrcs_dqlr.sh

echo "RSC ERASER, SWAP-------------------------------------------"

./scripts/probius/rsc_memory_eraser_swaps.sh

echo "RSC ERASER, DQLR-------------------------------------------"

./scripts/probius/rsc_memory_eraser_dqlr.sh
