#!/bin/sh

CODE=$1
DECODER=$2
JOBNAME=$(basename $CODE)

echo "#!/bin/bash
#SBATCH -J${JOBNAME}_MZ_5e-4_1e-3
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N64 --ntasks-per-node=8
#SBATCH --mem-per-cpu=2G
#SBATCH -t12:00:00
#SBATCH -qinferno
#SBATCH -oscripts/protean/evals/slurm/logs/%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

cd \$SLURM_SUBMIT_DIR
./scripts/protean/evals/compute_ber_5e-4_1e-3.sh ${CODE} ${DECODER} -mz" > "scripts/protean/evals/slurm/r${JOBNAME}_5e-4_1e-3_mz.sbatch"

echo "#!/bin/bash
#SBATCH -J${JOBNAME}_MX_5e-4_1e-3
#SBATCH --account=gts-mqureshi4-rg
#SBATCH -N64 --ntasks-per-node=8
#SBATCH --mem-per-cpu=2G
#SBATCH -t12:00:00
#SBATCH -qinferno
#SBATCH -oprotean/evals/slurm/logs/%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=suhaskvittal@gatech.edu

cd \$SLURM_SUBMIT_DIR
./scripts/protean/evals/compute_ber_5e-4_1e-3.sh ${CODE} ${DECODER} -mx" > "scripts/protean/evals/slurm/r${JOBNAME}_5e-4_1e-3_mx.sbatch"
