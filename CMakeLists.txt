# Declare CMAKE version 
cmake_minimum_required(VERSION 3.20.2)
project(QontraSim VERSION 0.1)
# Declare CPP Version
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(FLEX REQUIRED)
find_package(BISON REQUIRED)

BISON_TARGET(ASMParser ${CMAKE_SOURCE_DIR}/include/parsing/asm/parser.y 
    ${CMAKE_SOURCE_DIR}/include/parsing/asm/parser.tab.c
    DEFINES_FILE ${CMAKE_SOURCE_DIR}/include/parsing/asm/parser.tab.h)
FLEX_TARGET(ASMLexer ${CMAKE_SOURCE_DIR}/include/parsing/asm/lexer.lex 
    ${CMAKE_SOURCE_DIR}/include/parsing/asm/lex.yy.c)
ADD_FLEX_BISON_DEPENDENCY(ASMLexer ASMParser)

BISON_TARGET(SDLParser ${CMAKE_SOURCE_DIR}/include/parsing/sdl/parser.y
    ${CMAKE_SOURCE_DIR}/include/parsing/sdl/parser.tab.c
    DEFINES_FILE ${CMAKE_SOURCE_DIR}/include/parsing/sdl/parser.tab.h)
FLEX_TARGET(SDLLexer ${CMAKE_SOURCE_DIR}/include/parsing/sdl/lexer.lex 
    ${CMAKE_SOURCE_DIR}/include/parsing/sdl/lex.yy.c)
ADD_FLEX_BISON_DEPENDENCY(SDLLexer SDLParser)

set(QONTRA_FILES
    # Top-Level
    src/experiments.cpp
    src/instruction.cpp
    # Decoders
    src/decoder/mwpm.cpp
    src/decoder/window.cpp
    # Parsing
    src/parsing/cmd.cpp
    src/parsing/asm.cpp
    ${BISON_ASMParser_OUTPUTS}
    ${FLEX_ASMLexer_OUTPUTS}
#   src/parsing/sdl.cpp
#   ${BISON_SDLParser_OUTPUTS}
#   ${FLEX_SDLLexer_OUTPUTS}
    # Graphs
    src/graph/decoding_graph.cpp 
    src/graph/dependence_graph.cpp
#   src/graph/coupling_graph.cpp
    src/graph/lattice_graph.cpp
    src/graph/tanner_graph.cpp
    # Simulators
    src/sim/clifford_sim.cpp
    src/sim/frame_sim.cpp
    src/sim/manager.cpp
    src/sim/state_sim.cpp
    )

# Below are the remaining Qontra "extensions". The reason of I have
# not kept them with QONTRA_FILES is because not all systems
# are able to support them.
#
# For example, PACE does not support MLPACK natively.

if (COMPILE_NEURAL_DECODER)
    message("Will compile NEURAL DECODER.")
    find_package(Armadillo REQUIRED)
    set(QONTRA_FILES ${QONTRA_FILES} src/decoder/neural.cpp)
endif()

if (COMPILE_UNIVERSAL_SIM)
    message("Will compile UNIVERSAL SIM.")
    add_subdirectory(cytnx)
    set(QONTRA_FILES ${QONTRA_FILES} src/sim/universal_sim.cpp)
endif()

add_library(qontra ${QONTRA_FILES})

if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(qontra PUBLIC -Ofast -fno-strict-aliasing)
else()
    target_compile_options(qontra PUBLIC -ggdb3 -fno-strict-aliasing)
endif()

find_package(MPI REQUIRED)

target_include_directories(qontra PUBLIC "include" ${MPI_INCLUDE_PATH})

target_link_libraries(qontra PUBLIC libstim)
target_link_libraries(qontra PUBLIC libblossom5)
target_link_libraries(qontra PUBLIC ${MPI_CXX_LIBRARIES})

if (COMPILE_NEURAL_DECODER)
    target_include_directories(qontra PUBLIC ${ARMADILLO_INCLUDE_DIRS})
    target_link_libraries(qontra PUBLIC ${ARMADILLO_LIBRARIES})
endif()

if (COMPILE_UNIVERSAL_SIM)
    target_link_libraries(qontra PUBLIC cytnx)
endif()

add_subdirectory(stim)
add_subdirectory(blossom5)
add_subdirectory(pymatching)

add_executable(asm_analyzer src/main/asm_analyzer.cpp)
add_executable(converter src/main/converter.cpp)

target_link_libraries(asm_analyzer PRIVATE qontra)
target_link_libraries(converter PRIVATE qontra)

if (COMPILE_TESTS)
    message("Will compile TEST executables.")

    add_executable(clifford_sim_test src/main/tests/clifford_sim_test.cpp)
    add_executable(cytnx_contract_test src/main/tests/cytnx/contract_test.cpp)
    add_executable(tensor_sim_test src/main/tests/tensor_sim_test.cpp)

    target_link_libraries(clifford_sim_test PRIVATE qontra)
    target_link_libraries(cytnx_contract_test PRIVATE qontra)
    target_link_libraries(tensor_sim_test PRIVATE qontra)
endif()

if (COMPILE_MLDH)
    message("Will compile MLDH executables.")

    set(MLDH_FILES
        # Decoders
        src/mldh/block.cpp
        src/mldh/latency_sim.cpp)

    add_executable(mldh_memory src/main/personal/mldh/memory.cpp ${ASM_FILES})
    add_executable(mldh_wmemory src/main/personal/mldh/wmemory.cpp ${ASM_FILES})
    add_executable(mldh_bmemory src/main/personal/mldh/bmemory.cpp ${ASM_FILES} ${MLDH_FILES})

    add_executable(mldh_timing src/main/personal/mldh/timing.cpp ${MLDH_FILES})

    target_link_libraries(mldh_memory PRIVATE qontra)
    target_link_libraries(mldh_wmemory PRIVATE qontra)
    target_link_libraries(mldh_bmemory PRIVATE qontra)
    target_link_libraries(mldh_timing PRIVATE qontra)

    target_link_libraries(mldh_bmemory PRIVATE libpym)
endif()
