#!/bin/sh

CODE=$1
DECODER=$2
MEMORY=$3

folder="../data/protean/${CODE}/v4.4"
pmin=5e-4
pmax=1e-3

cd Release
#module load gcc/12.1.0-qgxpzk
#module load mvapich2/2.3.7-733lcv
echo "Running ${DECODER} for ${folder}"
mkdir "${folder}/output"
echo "mkdir ${folder}/output"
mpirun -np 144 pr_base_memory "${folder}" --decoder "${DECODER}" --e 50 "${MEMORY}" --pmin 5e-4 --pmax 1e-3 -fix-error

