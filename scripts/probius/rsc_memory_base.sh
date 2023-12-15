#!/bin/sh

trace_folder="../data/probius/surface_code_memory/${POLICY_NAME}/traces"
model_folder="../data/probius/surface_code_memory/${POLICY_NAME}/stim"
output_file="../data/probius/${POLICY_NAME}.csv"

cd build

rm -rf $trace_folder
rm -rf $model_folder

mkdir -p $model_folder

# Make traces
for d in $DISTANCES
do
    rounds=$(( d * ROUND_MULTIPLIER ))
    echo "tracing d = ${d}, rounds=${rounds}"
    for p in $ERRORS
    do
        name="memory_z_d${d}_p${p}"
        mkdir -p "${trace_folder}/${name}"
        echo "    p = ${p}"
        echo "    ---> tracing......."

        mpirun -np $PROC ./probius_tracer --d $d \
                                           --p $p \
                                           --shots $SHOTS \
                                           --rounds $rounds \
                                           --stim-out "${model_folder}/${name}.stim" \
                                           --trace-out "${trace_folder}/${name}" \
                                           --stats-out "${STATS_FOLDER}/${name}.csv" \
                                           $LEAKAGE_OPTIONS

        echo "    ---> decoding......."
        mpirun -np $PROC ./probius_memory --error-model "${model_folder}/${name}.stim" \
                                            --traces "${trace_folder}/${name}" \
                                            --out $output_file
    done
done

cd ..
