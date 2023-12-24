#!/bin/sh
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -Jrsc_mem_small_d
#SBATCH -N16 --ntasks-per-node=8
#SBATCH --mem-per-cpu=2G
#SBATCH -t2:00:00
#SBATCH -qinferno
#SBATCH -oslurm/probius/logs/rsc_mem_small_d.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

cd $SLURM_SUBMIT_DIR    # Should be repository root directory.

(
    export DISTANCES="3 5 7"
    export ERRORS="1e-4 5e-4 1e-3"
    export ROUND_MULTIPLIER=10

    export SYSTEM="pace"

    source scripts/probius_all.sh 128 100000
)
