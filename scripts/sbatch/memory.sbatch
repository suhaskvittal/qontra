#!/bin/bash

#SBATCH -Jmemory
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N8 --ntasks-per-node=16
#SBATCH --mem-per-cpu=4GB
#SBATCH -t6:00:00 
#SBATCH -qinferno
#SBATCH -osbatch/logs/memory.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

module load openmpi

cd build
for d in 7 9 11 13 15 17 19 21
do
    echo "Running bmemory for d = ${d}"
    for p in 1e-3
    do
        blk=$(( d/2 - 1 ))
        echo "\tParams: p = ${p}, r = ${d}, blk = ${blk}"
        srun ./mldh_bmemory --d ${d} --p ${p} --r ${d} --blk ${blk} --out ../data/pym_mem.csv --shots 1000000\
                            --timing-outdir ../data/timing/pym/d${d}_p${p} -pym
        srun ./mldh_bmemory --d ${d} --p ${p} --r ${d} --blk ${blk} --out ../data/astrea_pym_mem.csv --shots 1000000\
                            --timing-outdir ../data/timing/astrea_pym/d${d}_p${p}
    done
done
