# Author:   Suhas Vittal
# date:     25 December 2023

set(QONTRA_FILES
    # Top-Level
    src/qontra/experiments.cpp
    src/qontra/isa.cpp
    src/qontra/tables.cpp
    # Extensions
    src/qontra/ext/qes.cpp
    src/qontra/ext/stim.cpp
    # Decoders
    src/qontra/decoder/mwpm.cpp
    # Graphs
    src/qontra/graph/decoding_graph.cpp 
    src/qontra/graph/tanner_graph.cpp
    # Simulators
    src/qontra/sim/base/frame_sim.cpp
    src/qontra/sim/base/state_sim.cpp
    src/qontra/sim/full_system_sim/errors.cpp
    src/qontra/sim/full_system_sim/instructions.cpp
    src/qontra/sim/full_system_sim/run.cpp
    )

# Each extension may have its own source files. So, we will need to update
# the source files based on that.

if (COMPILE_NEURAL_DECODER)
    message(STATUS "Will compile neural network decoder. "
                    "Note that Armadillo and OpenMP must be available, "
                    "and MLPACK_INCLUDE_DIRS must be set.")
    find_package(Armadillo REQUIRED)
    find_package(OpenMP REQUIRED)
    set(QONTRA_FILES ${QONTRA_FILES} src/qontra/decoder/neural.cpp)
endif()

find_package(MPI REQUIRED)

add_library(qontra ${QONTRA_FILES})
target_compile_options(qontra PRIVATE ${COMPILE_OPTIONS} -fPIC)

if (L1D_CACHE_LINE_SIZE)
    target_compile_definitions(qontra PUBLIC L1D_CACHE_LINE_SIZE=${L1D_CACHE_LINE_SIZE})
endif()

if (COMPILE_PYMATCHING)
    target_compile_definitions(qontra PUBLIC QONTRA_PYMATCHING_ENABLED) 
endif()

if (COMPILE_CHROMOBIUS)
    target_compile_definitions(qontra PUBLIC QONTRA_CHROMOBIUS_ENABLED)
endif()
