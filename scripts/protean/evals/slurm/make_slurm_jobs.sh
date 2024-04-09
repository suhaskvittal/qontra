#!/bin/sh

CODE=$1
DECODER=$2

echo "#!/bin/sh

      #SBATCH -J${CODE}_MZ
      #SBATCH --account=gts-mqureshi4-rg
      #SBATCH -N64 --ntasks-per-node=8
      #SBATCH --mem-per-cpu=2G
      #SBATCH -t12:00:00
      #SBATCH -qinferno
      #SBATCH -oprotean/evals/slurm/logs/${CODE}_mz.out
      #SBATCH --mail-type=BEGIN,END,FAIL
      #SBATCH --mail-user=suhaskvittal@gatech.edu
      
      cd \$SLURM_SUBMIT_DIR
      ./scripts/protean/evals/compute_ber_5e-4_1e-3.sh ${CODE} ${DECODER} -mz" > "scripts/protean/evals/slurm/r${CODE}_mz.sbatch"

echo "#!/bin/sh

      #SBATCH -J${CODE}_MX
      #SBATCH --account=gts-mqureshi4-rg
      #SBATCH -N64 --ntasks-per-node=8
      #SBATCH --mem-per-cpu=2G
      #SBATCH -t12:00:00
      #SBATCH -qinferno
      #SBATCH -oprotean/evals/slurm/logs/${CODE}_mx.out
      #SBATCH --mail-type=BEGIN,END,FAIL
      #SBATCH --mail-user=suhaskvittal@gatech.edu
      
      cd \$SLURM_SUBMIT_DIR
      ./scripts/protean/evals/compute_ber_5e-4_1e-3.sh ${CODE} ${DECODER} -mx" > "scripts/protean/evals/slurm/r${CODE}_mx.sbatch"

