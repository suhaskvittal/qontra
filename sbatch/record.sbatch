#!/bin/bash

#SBATCH -Jrecord
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N8 --ntasks-per-node=4
#SBATCH --mem-per-cpu=4GB
#SBATCH -t00:30:00
#SBATCH -qinferno
#SBATCH -osbatch/logs/record.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

module load openmpi

cd build
d=11
blk=4
p=1e-3

srun ./mldh_bmemory --d ${d} --p ${p} --r ${d} --blk ${blk} --out ../data/bmemory.csv --shots 10000000\
                    --timing-outdir ../data/timing/d${d}_p${p} --syndrome-outdir ../data/syndromes/d${d}_p${p}\
                    --syndrome-record-hw 0
