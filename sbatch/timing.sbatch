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
for d in 7 9 11 13 15 17 19 21
do
    echo "Running timing for d = ${d}"
    for p in 1e-3
    do
        echo "\tParams: p = ${p}, r = ${d}, blk = ${blk}"
        t=$(( 1000*d ))
        srun ./mldh_timing --in ../data/timing/pym/d${d}_p${p} --out ../data/pym_timing.csv --t $t
#       srun ./mldh_timing --in ../data/timing/astrea_pym/d${d}_p${p} --out ../data/astrea_pym_timing.csv --t $t
    done
done
