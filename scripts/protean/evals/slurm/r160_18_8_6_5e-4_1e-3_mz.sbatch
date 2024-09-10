#!/bin/bash
#SBATCH -J160_18_8_6_MZ_5e-4_1e-3
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N64 --ntasks-per-node=8
#SBATCH --mem-per-cpu=2G
#SBATCH -t12:00:00
#SBATCH -qinferno
#SBATCH -oscripts/protean/evals/slurm/logs/%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

cd $SLURM_SUBMIT_DIR
./scripts/protean/evals/compute_ber_5e-4_1e-3.sh hysc/4_5/160_18_8_6 mwpm -mz
