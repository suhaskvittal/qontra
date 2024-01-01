# Author:   Suhas Vittal
# date:     25 December 2023

find_package(FLEX REQUIRED)
find_package(BISON REQUIRED)

# Setup Bison and Flex files.

BISON_TARGET(ASMParser ${CMAKE_SOURCE_DIR}/include/parsing/asm/parser.y 
    ${CMAKE_SOURCE_DIR}/include/parsing/asm/parser.tab.c
    DEFINES_FILE ${CMAKE_SOURCE_DIR}/include/parsing/asm/parser.tab.h)
FLEX_TARGET(ASMLexer ${CMAKE_SOURCE_DIR}/include/parsing/asm/lexer.lex 
    ${CMAKE_SOURCE_DIR}/include/parsing/asm/lex.yy.c)
ADD_FLEX_BISON_DEPENDENCY(ASMLexer ASMParser)

set(QONTRA_FILES
    # Top-Level
    src/experiments.cpp
    src/instruction.cpp
    src/tables.cpp
    # Utility and Definitions
    src/defs/filesystem.cpp
    # Decoders
    src/decoder/mwpm.cpp
    # Parsing
    src/parsing/cmd.cpp 
    src/parsing/asm.cpp
    ${BISON_ASMParser_OUTPUTS}
    ${FLEX_ASMLexer_OUTPUTS}
    # Graphs
    src/graph/decoding_graph.cpp 
    src/graph/lattice_graph.cpp
    src/graph/tanner_graph.cpp
    # Simulators
    src/sim/enumerator.cpp
    src/sim/frame_sim.cpp
    src/sim/state_sim.cpp
    src/sim/memory_sim.cpp
    )

# Each extension may have its own source files. So, we will need to update
# the source files based on that.

if (COMPILE_MEMORY_SIM_EXT)
    message(STATUS "Will compile memory simulator extensions.")
    set(QONTRA_FILES
            ${QONTRA_FILES}
            src/sim/memory_sim_ext/lrc.cpp
            src/sim/memory_sim_ext/eraser.cpp)
endif()

if (COMPILE_NEURAL_DECODER)
    message(STATUS "Will compile neural network decoder. Note that Armadillo and OpenMP must be available.")
    find_package(Armadillo REQUIRED)
    find_package(OpenMP REQUIRED)
    set(QONTRA_FILES ${QONTRA_FILES} src/decoder/neural.cpp)
endif()

find_package(MPI REQUIRED)

add_library(qontra ${QONTRA_FILES})
target_compile_options(qontra PRIVATE ${COMPILE_OPTIONS} -fPIC)

if (L1D_CACHE_LINE_SIZE)
    target_compile_definitions(qontra PUBLIC L1D_CACHE_LINE_SIZE=${L1D_CACHE_LINE_SIZE})
endif()

if (COMPILE_MEMORY_SIM_EXT)
    target_compile_definitions(qontra PUBLIC QONTRA_MEMORY_SIM_EXT_ENABLED)
endif()

if (COMPILE_PYMATCHING)
    target_compile_definitions(qontra PUBLIC QONTRA_PYMATCHING_ENABLED) 
endif()

if (COMPILE_CHROMOBIUS)
    target_compile_definitions(qontra PUBLIC QONTRA_CHROMOBIUS_ENABLED)
endif()
