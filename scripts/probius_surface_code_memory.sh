#!/bin/sh

proc=$1
shots=$2

trace_folder=../data/probius/surface_code_memory/traces
model_folder=../data/probius/surface_code_memory/stim
output_file=../data/probius/surface_code_memory.csv

cd build
make -j8

rm -rf $trace_folder
rm -rf $model_folder
mkdir -p $model_folder

# Make traces
for d in 3 5 7 9
do
    rounds=$((d+1))
    echo "d = ${d}, rounds=${rounds}"
    for p in 8e-4 9e-4 1e-3 2e-3 3e-3
    do
        name="memory_z_d${d}_p${p}"
        mkdir -p "${trace_folder}/${name}"
        echo "    p = ${p} | stim = ${model_folder}/${name}, trace = ${trace_folder}/${name}"
        mpirun -np $proc ./probius_tracer --d $d \
                                           --p $p \
                                           --shots $shots \
                                           --rounds $rounds \
                                           --stim-out "${model_folder}/${name}.stim" \
                                           --trace-out "${trace_folder}/${name}"
    done
done

# Run decoder
for d in 3 5 7 9
do
    for p in 8e-4 9e-4 1e-3 2e-3 3e-3
    do
        name="memory_z_d${d}_p${p}"
        mpirun -np $proc ./probius_memory --error-model "${model_folder}/${name}.stim" \
                                            --traces "${trace_folder}/${name}" \
                                            --out $output_file
    done
done
