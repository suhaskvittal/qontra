#!/bin/sh

PROC=$1
SHOTS=$2

# Setup environmental variables here.
# 
# First, we will identify the system we are using (user specified).
# Depending on that system's quirks, we need to modify our setup
# accordingly.

if [ -z "${SYSTEM}" ]; then
    MPI_CMD="mpirun -np ${PROC}"
    TRACER_BUILD_DIRECTORY=build
    MEMORY_BUILD_DIRECTORY=build
else
    if [ "${SYSTEM}" == "pace" ]; then
        echo "Setting up config for PACE..."

        MPI_CMD="srun"
        TRACER_BUILD_DIRECTORY=tracer_build
        MEMORY_BUILD_DIRECTORY=memory_build
        # Load any necessary modules here.
        module load openmpi
    fi
fi

export MPI_CMD
export TRACER_BUILD_DIRECTORY
export MEMORY_BUILD_DIRECTORY

# Now, we will determine if we are tracing, decoding, or both.

if [ -z "${MODE}" ]; then
    TRACE="yes"
    DECODE="yes"
elif [ "${MODE}" == "trace" ]; then
    TRACE="yes"
    DECODE="no"
else
    TRACE="no"
    DECODE="yes"
fi

export TRACE
export DECODE

echo "configuration: tracing = ${TRACE}, decoding = ${DECODE} --------------"

# Finally, export any other environmental variables we will need.

export SHOTS

if [ -z "${DISTANCES}" ]; then
    export DISTANCES="3 5 7"
fi

if [ -z "${ERRORS}" ]; then
    export ERRORS="1e-3 8e-4 6e-4"
fi

if [ -z "${ROUND_MULTIPLIER}" ]; then
    export ROUND_MULTIPLIER=5
fi

export STATS_FOLDER="../data/probius/surface_code_memory/simstats"

cd $TRACER_BUILD_DIRECTORY
mkdir -p $STATS_FOLDER
cd ..

