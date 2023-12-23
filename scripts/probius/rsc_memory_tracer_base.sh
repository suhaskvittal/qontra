#!/bin/sh

trace_folder="../data/probius/surface_code_memory/${POLICY_NAME}/traces"
model_folder="../data/probius/surface_code_memory/${POLICY_NAME}/stim"

cd $TRACER_BUILD_DIRECTORY

rm -rf $trace_folder
rm -rf $model_folder

mkdir -p $model_folder

for d in $DISTANCES
do
    rounds=$(( d * ROUND_MULTIPLIER ))
    echo "tracing d = ${d}, rounds = ${rounds}"
    for p in $ERRORS 
    do
        name="memory_z_d${d}_p${p}"
        mkdir -p "${trace_folder}/${name}"
        echo "    p = ${p}"

        $MPI_CMD ./probius_tracer --d $d \
                                   --p $p \
                                   --shots $SHOTS \
                                   --rounds $rounds \
                                   --stim-out "${model_folder}/${name}.stim" \
                                   --trace-out "${trace_folder}/${name}" \
                                   --stats-out "${STATS_FOLDER}/${name}.csv" \
                                   $LEAKAGE_OPTIONS
    done
done
