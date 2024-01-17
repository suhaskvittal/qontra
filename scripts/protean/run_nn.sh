#!/bin/sh

CODE=$1

: ${SHOTS=1000000}
: ${MPI_CMD=mpirun -np 8}
: ${MODEL_FILE=nn.bin}
: ${MODEL_FLAG=-nn}

cd Release

ARCH_FOLDER_PREFIX=../data/protean/${CODE}

for version in "1" "2.1" "2.2" "3.1" "3.2"
do
    echo "version ${version}"
    for p in 1e-5 2.5e-5 5e-5 7.5e-5 1e-4 2.5e-4 5e-4 7.5e-4 1e-3
    do
        echo "   p = ${p}"
        schedule_file=${ARCH_FOLDER_PREFIX}/v${version}/ext.qes
        output_file=${ARCH_FOLDER_PREFIX}/v${version}/memory.txt
        model_file=${ARCH_FOLDER_PREFIX}/v${version}/${MODEL_FILE}
        $MPI_CMD ./pr_nn_memory --qes $schedule_file \
                                --out $output_file \
                                --model $model_file \
                                --s $SHOTS \
                                --p $p \
                                ${MODEL_FLAG}
    done
done
