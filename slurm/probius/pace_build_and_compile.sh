#!/bin/sh
#
# INTENDED FOR USE WITH GEORGIA TECH'S PACE COMPUTING CLUSTER.
#
# Run from base directory of repository.

# Load modules

module load openmpi
module load cmake

total_proc=$1

rm -rf tracer_build
rm -rf memory_build

mkdir tracer_build
mkdir memory_build

# Make tracer executable first.

echo "TRACER -------------------------";

cd tracer_build
echo "---> building...";
cmake .. -DCMAKE_BUILD_TYPE=Release \
            -DCOMPILE_PROBIUS=ON \
            -DCOMPILE_MEMORY_SIM_EXT=ON \
            -DLINK_LEMON=ON \
            -DLEMON_INSTALL_DIR=/storage/coda1/p-mqureshi4/0/svittal8/pkgs \
            -DCOMPILE_PROBIUS_TRACER_ONLY=ON
echo "---> compiling...";
make -j${total_proc}

# Make memory executable now.

echo "MEMORY -------------------------";

cd ../memory_build
echo "---> building...";
cmake .. -DCMAKE_BUILD_TYPE=Release \
            -DCOMPILE_PROBIUS=ON \
            -DCOMPILE_PYMATCHING=ON \
            -DCOMPILE_PROBIUS_MEMORY_ONLY=ON
echo "---> compiling...";
make -j${total_proc}

