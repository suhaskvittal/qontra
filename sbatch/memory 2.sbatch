#!/bin/bash

#SBATCH -Jmemory
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N8 --ntasks-per-node=16
#SBATCH --mem-per-cpu=4GB
#SBATCH -t00:30:00
#SBATCH -qinferno
#SBATCH -osbatch/logs/memory.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

module load openmpi

cd build
for d in 5 7 9 11
do
    echo "Running bmemory for d = ${d}"
    for p in 3e-3 2e-3 1e-3 9e-4 8e-4 7e-4
    do
        blk=$(( d/2 - 1 ))
        echo "\tParams: p = ${p}, r = ${d}, blk = ${blk}"
        srun ./mldh_bmemory --d ${d} --p ${p} --r ${d} --blk ${blk} --out ../data/bmemory.csv --shots 10000000\
                            --timing-outdir ../data/timing/d${d}_p${p}
    done
done
