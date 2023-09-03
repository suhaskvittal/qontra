#!/bin/bash

#SBATCH -Jtiming
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N4 --ntasks-per-node=4
#SBATCH --mem-per-cpu=4GB
#SBATCH -t00:30:00
#SBATCH -qinferno
#SBATCH -osbatch/logs/timing.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

module load openmpi

cd build
for d in 5 7 9 11
do
    echo "Running timing for d = ${d}"
    for p in 3e-3 2e-3 1e-3 9e-4 8e-4 7e-4
    do
        echo "\tParams: p = ${p}, r = ${d}, blk = ${blk}"
        srun ./mldh_timing --in ../data/timing/d${d}_p${p} --out ../data/timing.csv
    done
done
