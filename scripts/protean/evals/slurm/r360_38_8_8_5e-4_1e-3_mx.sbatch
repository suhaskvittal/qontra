#!/bin/bash
#SBATCH -J360_38_8_8_MX_5e-4_1e-3
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N64 --ntasks-per-node=8
#SBATCH --mem-per-cpu=12G
#SBATCH -t24:00:00
#SBATCH -qinferno
#SBATCH -oprotean/evals/slurm/logs/%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

cd $SLURM_SUBMIT_DIR
./scripts/protean/evals/compute_ber_5e-4_1e-3.sh hysc/4_5/360_38_8_8 mwpm -mx
