#!/bin/bash
#SBATCH -J480_148_6_6_MZ_5e-4_1e-3
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N64 --ntasks-per-node=8
#SBATCH --mem-per-cpu=8G
#SBATCH -t24:00:00
#SBATCH -qinferno
#SBATCH -oscripts/protean/evals/slurm/logs/%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

cd $SLURM_SUBMIT_DIR
./scripts/protean/evals/compute_ber_5e-4_1e-3.sh hycc/4_10/480_148_6_6 restriction -mz
