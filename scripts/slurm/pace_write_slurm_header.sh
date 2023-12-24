#!/bin/sh
#
# INTENDED FOR USE WITH GEORGIA TECH's PACE COMPUTING CLUSTER

job_name=$1

nodes=$2
cores=$3    # Per node
memory=$4   # Per core
walltime=$5

# Setup file variables here.

output_file="${job_name}.sbatch"
log_file="logs/${job_name}.log"

mkdir logs

echo "#!/bin/sh" >> $output_file
echo "#SBATCH --account=gts-mqureshi4-rg" >> $output_file

echo "#SBATCH -J${job_name}" >> $output_file

echo "#SBATCH -N${nodes} --ntasks-per-node=${cores}" >> $output_file
echo "#SBATCH --mem-per-cpu=${memory}" >> $output_file
echo "#SBATCH -t${walltime}" >> $output_file

echo "#SBATCH -qinferno" >> $output_file
echo "#SBATCH -o${log_file}" >> $output_file
echo "#SBATCH --mail-type=BEGIN,END,FAIL" >> $output_file
echo "#SBATCH --mail-user=suhaskvittal@gatech.edu" >> $output_file

echo "" >> $output_file
echo "cd \$SLURM_SUBMIT_DIR" >> $output_file
