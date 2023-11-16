#!/bin/bash

#SBATCH -Jcompile_mldh
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N1 --ntasks-per-node=8
#SBATCH --mem-per-cpu=256M
#SBATCH -t00:05:00
#SBATCH -qinferno
#SBATCH -osbatch/logs/compile.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

module load cmake
module load openmpi
module load gcc

rm -rf build
mkdir build && cd build

cmake .. -DCMAKE_BUILD_TYPE=Release -DCOMPILE_MLDH=On
make -j8
