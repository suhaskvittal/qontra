#!/bin/sh

cd build

shots=100000000

for p in 1e-5 2e-5 3e-5
#for p in 4e-5 5e-5 6e-5 7e-5 8e-5 9e-5
#for p in 1e-4 2e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 1e-3
do
    echo "running for p = ${p}"
    mpirun -np 6 ./memory --asm ../projects/qontra/asm/qec/hycc/mzd4_both_at_once_no_flags.asm --out ../projects/hycc.csv --p $p --shots $shots --model-file ../projects/models/hycc_both_at_once_no_flags.bin
#   mpirun -np 6 ./memory --asm ../projects/qontra/asm/qec/hycc/mzd4_both_at_once.asm --out ../projects/hycc.csv --p $p --shots $shots --model-file ../projects/models/hycc_both_at_once.bin
#   mpirun -np 6 ./memory --asm ../projects/qontra/asm/qec/hycc/mzd4_generated_v1.asm --out ../projects/hycc.csv --p $p --shots $shots --model-file ../projects/models/hycc_v1.bin
#   mpirun -np 6 ./memory --asm ../projects/qontra/asm/qec/hycc/mzd4_generated_v2.asm --out ../projects/hycc.csv --p $p --shots $shots --model-file ../projects/models/hycc_v2.bin
#    mpirun -np 6 ./memory --asm ../projects/qontra/asm/qec/hycc/mzd4_generated_v3.asm --out ../projects/hycc.csv --p $p --shots $shots --model-file ../projects/models/hycc_v3.bin
#   mpirun -np 6 ./memory --asm ../projects/qontra/asm/qec/hycc/mzd4_generated_v5.asm --out ../projects/hycc.csv --p $p --shots $shots --model-file ../projects/models/hycc_v5.bin
done
