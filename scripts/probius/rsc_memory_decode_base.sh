#!/bin/sh

trace_folder="../data/probius/surface_code_memory/${POLICY_NAME}/traces"
model_folder="../data/probius/surface_code_memory/${POLICY_NAME}/stim"
output_file="../data/probius/${POLICY_NAME}.csv"

cd $MEMORY_BUILD_DIRECTORY

for d in $DISTANCES
do
    rounds=$(( d * ROUND_MULTIPLIER ))
    echo "decoding d = ${d}, rounds = ${rounds}"
    for p in $ERRORS
    do
        name="memory_z_d${d}_p${p}"
        echo "    p = ${p}"

        $MPI_CMD ./probius_memory --error-model "${model_folder}/${name}.stim" \
                                    --traces "${trace_folder}/${name}" \
                                    --out $output_file
    done
done
