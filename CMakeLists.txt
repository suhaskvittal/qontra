cmake_minimum_required(VERSION 3.20.2)
project(QontraSim VERSION 0.1)

#
# QontraSim CMAKE Options:
#   (1) COMPILE_MEMORY_SIM_EXT : compiles memory simulator extensions (i.e. LRCs).
#   (2) COMPILE_NEURAL_DECODER : compiles the neural network decoder written with MLPACK.
#   (3) COMPILE_PYMATCHING : compiles PyMatching v2.
#   (4) LINK_LEMON : user must specify one of:
#                       (i) LEMON_LIBRARY_DIR and LEMON_INCLUDE_DIR
#                       (ii) LEMON_INSTALL_PATH, which sets 
#                               LEMON_LIBRARY_DIR=${LEMON_INSTALL_PATH}/lemon/lib
#                               LEMON_INCLUDE_DIR=${LEMON_INSTALL_PATH}/lemon/include
#                    The user must compile LEMON on their end.
#

# Declare CPP Version
set(CMAKE_CXX_STANDARD 17)
# Some libraries may require C++ versions higher than 17.
#   -- Pymatching uses <bit>, which requires C++20.
if (COMPILE_PYMATCHING)
    set(CMAKE_CXX_STANDARD 20)
endif()

if (COMPILE_CHROMOBIUS)
    set(CMAKE_CXX_STANDARD 20)
    set(COMPILE_PYMATCHING On)  # PyMatching is a dependency.
endif()

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
#   src/decoder/restriction.cpp
    # Parsing
    src/parsing/cmd.cpp src/parsing/asm.cpp
    ${BISON_ASMParser_OUTPUTS}
    ${FLEX_ASMLexer_OUTPUTS}
    # Graphs
    src/graph/decoding_graph.cpp 
    src/graph/lattice_graph.cpp
    src/graph/tanner_graph.cpp
    # Simulators
#   src/sim/clifford_sim.cpp        # DISABLED FOR NOW -- RE-ENABLE LATER
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

if (COMPILE_MEMORY_SIM_EXT)
    message(STATUS "Will compile memory simulator extensions.")
    set(QONTRA_FILES
            ${QONTRA_FILES}
            src/sim/memory_sim_ext/lrc.cpp
            src/sim/memory_sim_ext/eraser.cpp)
endif()

if (COMPILE_PYMATCHING)
    message(STATUS "Will compile PyMatching.")
    add_subdirectory(dependencies/pymatching)
endif()

if (COMPILE_CHROMOBIUS)
    message(STATUS "Will compile Chromobius.")
    add_subdirectory(dependencies/chromobius)
endif()

if (COMPILE_NEURAL_DECODER)
    message(STATUS "Will compile NEURAL DECODER.")
    find_package(Armadillo REQUIRED)
    set(QONTRA_FILES ${QONTRA_FILES} src/decoder/neural.cpp)
endif()

find_package(MPI REQUIRED)

add_library(qontra ${QONTRA_FILES})

if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(qontra PUBLIC -O3 -fno-strict-aliasing)
else()
    target_compile_options(qontra PUBLIC -ggdb3 -fno-strict-aliasing)
endif()

if (COMPILE_MEMORY_SIM_EXT)
    target_compile_definitions(qontra PUBLIC QONTRA_MEMORY_SIM_EXT_ENABLED)
endif()

target_include_directories(qontra PUBLIC "include" ${MPI_INCLUDE_PATH} "/usr/local/include")

target_link_libraries(qontra PUBLIC libstim)
target_link_libraries(qontra PUBLIC libblossom5)
target_link_libraries(qontra PUBLIC ${MPI_CXX_LIBRARIES})

if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(qontra PUBLIC "stdc++fs")
endif()

if (COMPILE_PYMATCHING)
    target_link_libraries(qontra PUBLIC libpym)
    target_compile_definitions(qontra PUBLIC QONTRA_PYMATCHING_ENABLED)
endif()

if (COMPILE_CHROMOBIUS)
    target_link_libraries(qontra PUBLIC libchromobius)
    target_compile_definitions(qontra PUBLIC QONTRA_CHROMOBIUS_ENABLED)
endif()

if (COMPILE_NEURAL_DECODER)
    message(STATUS "Adding Armadillo include and libs")
    target_include_directories(qontra PUBLIC ${ARMADILLO_INCLUDE_DIRS})
    target_link_libraries(qontra PUBLIC ${ARMADILLO_LIBRARIES})

    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        message(STATUS "Found OpenMP for the NN Decoder.")
        target_link_libraries(qontra PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()

if (LINK_LEMON)
    if (DEFINED LEMON_INSTALL_DIR)
        set(LEMON_LIBRARY_DIR "${LEMON_INSTALL_DIR}/lemon/lib")
        set(LEMON_INCLUDE_DIR "${LEMON_INSTALL_DIR}/lemon/include")
    endif()

    if (NOT DEFINED LEMON_LIBRARY_DIR)
        message(FATAL_ERROR "COMPILE_LEMON is set, but LEMON_LIBRARY_DIR is not defined.")
    elseif(NOT DEFINED LEMON_INCLUDE_DIR)
        message(FATAL_ERROR "COMPILE_LEMON is set, but LEMON_INCLUDE_DIR is not defined.")
    else()
        # Everything checks out. Link lemon.
        message(STATUS "Will link Lemon, LEMON_LIBRARY_DIR=${LEMON_LIBRARY_DIR}, LEMON_INCLUDE_DIR=${LEMON_INCLUDE_DIR}")
        add_library(lemon STATIC IMPORTED)
        set_target_properties(lemon PROPERTIES IMPORTED_LOCATION "${LEMON_LIBRARY_DIR}/libemon.a")
        target_link_libraries(qontra PUBLIC lemon)
        target_include_directories(qontra PUBLIC ${LEMON_INCLUDE_DIR})
    endif()
endif()

add_subdirectory(dependencies/stim)
add_subdirectory(dependencies/blossom5)

add_executable(asm_analyzer main/asm_analyzer.cpp)
add_executable(converter main/converter.cpp)
add_executable(memory main/memory.cpp)
add_executable(syndromes main/syndromes.cpp)

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

    add_executable(protean_scheduler_test main/protean/tests/scheduler.cpp)
    add_executable(protean_writing_test main/protean/tests/writing.cpp)
    target_link_libraries(protean_scheduler_test PRIVATE protean)
    target_link_libraries(protean_writing_test PRIVATE protean)
endif()

if (COMPILE_PROBIUS)
    # This is because PyMatching needs C++20, but LEMON does not work with
    # C++20 for whatever reason.
    if (COMPILE_PROBIUS_TRACER_ONLY)
        set(COMPILE_PROBIUS_TRACER ON)
        set(COMPILE_PROBIUS_MEMORY OFF)
    elseif (COMPILE_PROBIUS_MEMORY_ONLY)
        set(COMPILE_PROBIUS_TRACER OFF)
        set(COMPILE_PROBIUS_MEMORY ON)
    else()
        set(COMPILE_PROBIUS_TRACER ON)
        set(COMPILE_PROBIUS_MEMORY ON)
    endif()

    if (COMPILE_PROBIUS_TRACER)
        add_executable(probius_tracer main/probius/tracer.cpp)
        target_link_libraries(probius_tracer PRIVATE qontra)
    endif()

    if (COMPILE_PROBIUS_MEMORY)
        add_executable(probius_memory main/probius/memory.cpp)
        target_link_libraries(probius_memory PRIVATE qontra)
    endif()
endif()
