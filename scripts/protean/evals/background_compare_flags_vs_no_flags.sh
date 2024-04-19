#!/bin/sh

PROC=$1

CODE='hysc/4_5/60_8_6_4'

cd Release
mkdir -p "../data/protean/${CODE}/v1"
mkdir -p "../data/protean/${CODE}/v2"

./protean ../data/tanner/${CODE}.txt ../data/protean/${CODE}/v1 \
            --passes "rlb.rcr" \
            --s-rounds 4 \
            --max-conn 4 \
            -skip-render
./protean ../data/tanner/${CODE}.txt ../data/protean/${CODE}/v2 \
            --passes "fla.ral.rlb.rcr" \
            --s-rounds 4 \
            --max-conn 4 \
            -skip-render

for vx in 'v1' 'v2'
do
    mkdir ../data/protean/${CODE}/${vx}/output
    mpirun -np ${PROC} ./pr_base_memory ../data/protean/${CODE}/${vx} --decoder mwpm\
                --e 50 --pmin 5e-4 --pmax 1e-3 -fix-error
done
