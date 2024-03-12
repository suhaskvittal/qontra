#!/bin/sh

cd data

dt=`date '+%m_%d_%y_%H%M'`
fname="protean_${dt}.tar.gz"
echo "Sending data at date ${dt}"

tar -czvf $fname protean/
scp $fname $PACE:~/p-mqureshi4-0/qec/qontra/data/

