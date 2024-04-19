#!/bin/sh

PROC=$1

mkdir -p data/protean/rsc/output
mkdir -p data/protean/hexcc/output
cd Release

for d in 3 5 7
do
    for mem in "x" "z"
    do
        echo "Running d = ${d}, mem = ${mem}, rsc"
        mpirun -np ${PROC} ./pr_planar_memory ../data/protean/rsc/${mem}_d${d}.qes \
                    ../data/protean/rsc/output/${mem}_d${d}.csv \
                    --decoder mwpm --e 50 --pmin 5e-4 --pmax 1e-3

        echo "Running d = ${d}, mem = ${mem}, hexcc"
        mpirun -np ${PROC} ./pr_planar_memory ../data/protean/hexcc/${mem}_d${d}.qes \
                    ../data/protean/hexcc/output/${mem}_d${d}.csv \
                    --decoder restriction --e 50 --pmin 5e-4 --pmax 1e-3
    done
done
