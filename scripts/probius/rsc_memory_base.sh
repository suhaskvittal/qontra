#!/bin/sh

DISTANCES="3 5 7"
ERRORS="1e-3"
ROUND_MULTIPLIER=5

stats_folder="../data/probius/surface_code_memory/${POLICY_NAME}/simstats"
trace_folder="../data/probius/surface_code_memory/${POLICY_NAME}/traces"
model_folder="../data/probius/surface_code_memory/${POLICY_NAME}/stim"
output_file="../data/probius/${POLICY_NAME}.csv"

cd build

rm -rf $stats_folder
rm -rf $trace_folder
rm -rf $model_folder

mkdir -p $stats_folder
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
                                           --stats-out "${stats_folder}/${name}.csv" \
                                           $LEAKAGE_OPTIONS

        echo "    ---> decoding......."
        mpirun -np $PROC ./probius_memory --error-model "${model_folder}/${name}.stim" \
                                            --traces "${trace_folder}/${name}" \
                                            --out $output_file
    done
done

cd ..
