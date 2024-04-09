#!/bin/sh

CODE=$1
DECODER=$2
MEMORY=$3

folder="data/protean/$CODE/v4.4"
pmin=5e-4
pmax=1e-3

cd Release
srun ./pr_base_memory $folder --decoder $DECODER --pmin $pmin --pmax $pmax --e 50 -fix-error $MEMORY

