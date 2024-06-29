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
    src/qontra/decoder/concat_mwpm.cpp
    src/qontra/decoder/matching_base.cpp
    src/qontra/decoder/mobius.cpp
    src/qontra/decoder/mwpm.cpp
    src/qontra/decoder/restriction.cpp
    src/qontra/decoder/restriction/helpers.cpp
    src/qontra/decoder/restriction/incident_vertices.cpp
    # Graphs
    src/qontra/graph/decoding_graph.cpp 
    src/qontra/graph/decoding_graph/detectors.cpp
    src/qontra/graph/decoding_graph/distance.cpp
    src/qontra/graph/decoding_graph/edge_class.cpp
    src/qontra/graph/decoding_graph/helpers.cpp
    src/qontra/graph/decoding_graph/init.cpp
    src/qontra/graph/decoding_graph/rgb_only_lattice.cpp
    src/qontra/graph/decoding_graph/unified_lattice.cpp
    src/qontra/graph/tanner_graph.cpp
    src/qontra/graph/tanner_graph/code_distance.cpp
    src/qontra/graph/tanner_graph/io.cpp
    # Simulators
    src/qontra/sim/base/frame_sim.cpp
    src/qontra/sim/base/state_sim.cpp
    src/qontra/sim/base/clifford_sim.cpp
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
    set(QONTRA_FILES ${QONTRA_FILES} 
            src/qontra/decoder/neural.cpp)
endif()

find_package(MPI REQUIRED)

add_library(qontra ${QONTRA_FILES})
target_compile_options(qontra PRIVATE ${COMPILE_OPTIONS} -fPIC)

target_compile_definitions(qontra PUBLIC QONTRA_ISA_FILE="${QONTRA_ISA_FILE}")
#target_compile_definitions(qontra PUBLIC DECODER_PERF)

if (L1D_CACHE_LINE_SIZE)
    target_compile_definitions(qontra PUBLIC L1D_CACHE_LINE_SIZE=${L1D_CACHE_LINE_SIZE})
endif()

if (COMPILE_PYMATCHING)
    set(QONTRA_FILES ${QONTRA_FILES} src/qontra/decoder/pymatching.cpp)
    target_compile_definitions(qontra PUBLIC QONTRA_PYMATCHING_ENABLED) 
endif()

if (COMPILE_CHROMOBIUS)
    set(QONTRA_FILES ${QONTRA_FILES} src/qontra/decoder/chromobius.cpp)
    target_compile_definitions(qontra PUBLIC QONTRA_CHROMOBIUS_ENABLED)
endif()
