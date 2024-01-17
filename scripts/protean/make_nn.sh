#!/bin/sh

CODE=$1

: ${SHOTS=100000}
: ${EPOCHS=300}
: ${MODEL_FILE=nn.bin}
: ${MODEL_FLAG=-nn}

echo "training params: shots = ${SHOTS}, epochs = ${EPOCHS}"

cd Release

ARCH_FOLDER_PREFIX=../data/protean/${CODE}

for version in "1" "2.1" "2.2" "3.1" "3.2"
do
    schedule_file=${ARCH_FOLDER_PREFIX}/v${version}/ext.qes
    model_file=${ARCH_FOLDER_PREFIX}/v${version}/${MODEL_FILE}
    ./pr_nn_train --qes $schedule_file \
                    --out $model_file \
                    --s $SHOTS \
                    --e $EPOCHS \
                    --p 5e-4 \
                    ${MODEL_FLAG}
done

