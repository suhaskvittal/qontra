# Declare CMAKE version 
cmake_minimum_required(VERSION 3.20.2)
project(QontraSim VERSION 0.1)
# Declare CPP Version
set(CMAKE_CXX_STANDARD 17)
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
#   src/decoder/restriction.cpp
    # Parsing
    src/parsing/cmd.cpp
    src/parsing/asm.cpp
    ${BISON_ASMParser_OUTPUTS}
    ${FLEX_ASMLexer_OUTPUTS}
    # Graphs
    src/graph/decoding_graph.cpp 
    src/graph/dependence_graph.cpp
    src/graph/lattice_graph.cpp
    src/graph/tanner_graph.cpp
    # Simulators
    src/sim/clifford_sim.cpp
    src/sim/enumerator.cpp
    src/sim/frame_sim.cpp
    src/sim/state_sim.cpp
    src/sim/memory_sim.cpp
    )

# Compatability
if (APPLE)
    set(CMAKE_CXX_FLAGS "-no-pie")
endif()

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
    target_compile_options(qontra PUBLIC -O3 -fno-strict-aliasing)
else()
    target_compile_options(qontra PUBLIC -ggdb3 -fno-strict-aliasing)
endif()

find_package(MPI REQUIRED)

target_include_directories(qontra PUBLIC "include" ${MPI_INCLUDE_PATH} "/usr/local/include")

target_link_libraries(qontra PUBLIC libstim)
target_link_libraries(qontra PUBLIC libblossom5)
target_link_libraries(qontra PUBLIC ${MPI_CXX_LIBRARIES})

if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(qontra PUBLIC "stdc++fs")
endif()

if (COMPILE_NEURAL_DECODER)
    message("Adding Armadillo include and libs")
    target_include_directories(qontra PUBLIC ${ARMADILLO_INCLUDE_DIRS})
    target_link_libraries(qontra PUBLIC ${ARMADILLO_LIBRARIES})

    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        message("Found OpenMP for the NN Decoder.")
        target_link_libraries(qontra PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()

if (COMPILE_UNIVERSAL_SIM)
    target_link_libraries(qontra PUBLIC cytnx)
endif()

add_subdirectory(stim)
add_subdirectory(blossom5)

add_executable(asm_analyzer src/main/asm_analyzer.cpp)
add_executable(converter src/main/converter.cpp)
add_executable(memory src/main/memory.cpp)
add_executable(syndromes src/main/syndromes.cpp)

target_link_libraries(asm_analyzer PRIVATE qontra)
target_link_libraries(converter PRIVATE qontra)
target_link_libraries(memory PRIVATE qontra)
target_link_libraries(syndromes PRIVATE qontra)

if (COMPILE_PROTEAN)
    include(${CMAKE_SOURCE_DIR}/cmake/FindCPLEX.cmake)

    set(PROTEAN_FILES
        src/protean/scheduler.cpp
        src/protean/utils.cpp)
    add_library(protean ${PROTEAN_FILES})
    target_include_directories(protean PUBLIC ${CPLEX_INCLUDE_DIR})
    target_link_libraries(protean PUBLIC 
                            qontra 
                            ${CPLEX_LIB_PREFIX}/libcplex.a
                            ${CPLEX_LIB_PREFIX}/libilocplex.a)
    target_compile_options(protean PUBLIC -O3 -fno-strict-aliasing)

    add_executable(protean_scheduler_test src/main/protean/tests/scheduler.cpp)
    add_executable(protean_writing_test src/main/protean/tests/writing.cpp)
    target_link_libraries(protean_scheduler_test PRIVATE protean)
    target_link_libraries(protean_writing_test PRIVATE protean)
endif()

if (COMPILE_PROBIUS)
    add_executable(probius_tracer src/main/probius/tracer.cpp)
    add_executable(probius_memory src/main/probius/memory.cpp)
    target_link_libraries(probius_tracer PRIVATE qontra)
    target_link_libraries(probius_memory PRIVATE qontra)
endif()

if (COMPILE_TESTS)
    message("Will compile TEST executables.")

    add_executable(clifford_sim_test src/main/tests/clifford_sim_test.cpp)
    add_executable(cytnx_contract_test src/main/tests/cytnx/contract_test.cpp)
    add_executable(tensor_sim_test src/main/tests/tensor_sim_test.cpp)

    target_link_libraries(clifford_sim_test PRIVATE qontra)
    target_link_libraries(cytnx_contract_test PRIVATE qontra)
    target_link_libraries(tensor_sim_test PRIVATE qontra)
endif()
