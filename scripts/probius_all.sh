#!/bin/sh

proc=$1
shots=$2

(

export DISTANCES="3 5 7"
export ERRORS="1e-3 8e-4 6e-4"
export ROUND_MULTIPLIER=5
export STATS_FOLDER="../data/probius/surface_code_memory/simstats"

cd build
mkdir -p $STATS_FOLDER
cd ..

echo "RSC NO LEAKAGE---------------------------------------------"

#./scripts/probius/rsc_memory_no_leakage.sh $proc $shots

echo "RSC LEAKAGE, NO LRCS---------------------------------------"

#./scripts/probius/rsc_memory_leakage_no_lrcs.sh $proc $shots

echo "RSC ALWAYS LRCS, SWAP--------------------------------------"

./scripts/probius/rsc_memory_always_lrcs_swaps.sh $proc $shots

echo "RSC ALWAYS LRCS, DQLR--------------------------------------"

./scripts/probius/rsc_memory_always_lrcs_dqlr.sh $proc $shots

echo "RSC OPTIMAL LRCS, SWAP-------------------------------------"

./scripts/probius/rsc_memory_optimal_lrcs_swaps.sh $proc $shots

echo "RSC OPTIMAL LRCS, DQLR-------------------------------------"

./scripts/probius/rsc_memory_optimal_lrcs_dqlr.sh $proc $shots

echo "RSC ERASER, SWAP-------------------------------------------"

./scripts/probius/rsc_memory_eraser_swaps.sh $proc $shots

echo "RSC ERASER, DQLR-------------------------------------------"

./scripts/probius/rsc_memory_eraser_dqlr.sh $proc $shots

)
