#!/bin/sh

PROC=$1

CODE='hysc/4_5/60_8_6_4'

cd Release
mkdir -p "../data/protean/${CODE}"

for ps in 'rlb.rcr' 'fla.ral.rlb.rcr'
do
    ./protean ../data/tanner/{CODE}.txt ../data/protean/${CODE} \
                --passes $ps \
                --s-rounds 4 \
                --max-conn 4 \
                -skip-render
done
