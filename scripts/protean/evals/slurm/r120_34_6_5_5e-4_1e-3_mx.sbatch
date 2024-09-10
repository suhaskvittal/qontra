#!/bin/bash
#SBATCH -J120_34_6_5_MX_5e-4_1e-3
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N64 --ntasks-per-node=8
#SBATCH --mem-per-cpu=2G
#SBATCH -t12:00:00
#SBATCH -qinferno
#SBATCH -oprotean/evals/slurm/logs/%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

cd $SLURM_SUBMIT_DIR
./scripts/protean/evals/compute_ber_5e-4_1e-3.sh hysc/5_6/120_34_6_5 mwpm -mx
