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
    ${CMAKE_SOURCE_DIR}/include/parsing/asm/lex.yy.c
    DEFINES_FILE ${CMAKE_SOURCE_DIR}/include/parsing/asm/lex.tab.h)
ADD_FLEX_BISON_DEPENDENCY(ASMLexer ASMParser)

set(QONTRA_FILES
    # Top-Level
    src/experiments.cpp
    src/instruction.cpp
    # Decoders
    src/decoder/mwpm.cpp
    src/decoder/window.cpp
    # Parsing
    src/parsing/cmd.cpp
    # Graphs
    src/graph/decoding_graph.cpp 
    src/graph/coupling_graph.cpp
    src/graph/tanner_graph.cpp
    # Protean
    src/protean/compiler.cpp
    src/protean/proc3d.cpp
    # Simulation
    src/sim/clifford_sim.cpp
    src/sim/control_sim.cpp
    src/sim/frame_sim.cpp
    src/sim/state_sim.cpp
    )
# I've a had a lot of trouble with Bison/Flex files. Only
# include them if you need them.
set(ASM_FILES
    src/parsing/asm.cpp
    ${BISON_ASMParser_OUTPUTS}
    ${FLEX_ASMLexer_OUTPUTS}
    )

add_library(qontra ${QONTRA_FILES})
add_executable(memory src/main/memory_experiment.cpp)
add_executable(circuit src/main/build_circuit.cpp)
add_executable(stim_to_asm src/main/stim_to_asm.cpp ${ASM_FILES})

### Protean ###
add_executable(protean_main src/main/protean/main.cpp)

### TESTING ###
add_executable(ctrl src/main/test/ctrlsimtest.cpp ${ASM_FILES})
add_executable(asmtest src/main/test/asmparsetest.cpp ${ASM_FILES})
add_executable(perceptron src/main/test/perceptron_test.cpp)

if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(qontra PUBLIC -O3 -lemon -fno-strict-aliasing)
else()
    target_compile_options(qontra PUBLIC -ggdb3 -lemon -fno-strict-aliasing)
endif()

find_package(MPI REQUIRED)
find_package(LEMON REQUIRED)

# I don't know why I have to include /usr/local/include...
target_include_directories(qontra PUBLIC
    "include" 
    ${MPI_INCLUDE_PATH} 
    "/usr/local/include")

add_subdirectory(stim)
add_subdirectory(blossom5)

target_link_libraries(qontra PUBLIC libstim)
target_link_libraries(qontra PUBLIC libblossom5)
target_link_libraries(qontra PUBLIC ${MPI_CXX_LIBRARIES})

target_link_libraries(memory PRIVATE qontra)
target_link_libraries(circuit PRIVATE qontra)
target_link_libraries(stim_to_asm PRIVATE qontra)

target_link_libraries(ctrl PRIVATE qontra)
target_link_libraries(asmtest PRIVATE qontra)
target_link_libraries(perceptron PRIVATE qontra)

target_link_libraries(protean_main PRIVATE qontra)
