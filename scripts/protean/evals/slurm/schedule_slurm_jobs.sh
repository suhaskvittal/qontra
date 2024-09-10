#!/bin/sh
CODE=$1
sbatch "scripts/protean/evals/slurm/r${CODE}_5e-4_1e-3_mz.sbatch"
sbatch "scripts/protean/evals/slurm/r${CODE}_5e-4_1e-3_mx.sbatch"
