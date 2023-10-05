#!/bin/bash

#SBATCH -Jwindow
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N16 --ntasks-per-node=8
#SBATCH --mem-per-cpu=4GB
#SBATCH -t00:30:00
#SBATCH -qinferno
#SBATCH -osbatch/logs/window.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

module load openmpi

cd build

p=1e-3
shots=10000000

for d in 5 7 9 11
do
    echo "Running window for d = ${d}"
    rounds=$(( 3*d )) 

    comm1=1
    comm2=$(( d/2 ))
    comm3=$d

    for c in $comm1 $comm2 $comm3
    do
        srun ./mldh_wmemory --d ${d} --p ${p} --r ${rounds} --w ${d} --c ${c} --out ../data/window.csv --shots ${shots}
    done
done
