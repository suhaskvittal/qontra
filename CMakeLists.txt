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
    src/decoder/restriction.cpp
    src/decoder/neural.cpp
    # Parsing
    src/parsing/cmd.cpp
    # Graphs
    src/graph/decoding_graph.cpp 
    src/graph/dependence_graph.cpp
    src/graph/coupling_graph.cpp
    src/graph/lattice_graph.cpp
    src/graph/tanner_graph.cpp
    # Simulators
    src/sim/clifford_sim.cpp
    src/sim/frame_sim.cpp
    src/sim/manager.cpp
    src/sim/state_sim.cpp
    )
# I've a had a lot of trouble with Bison/Flex files. Only
# include them if you need them.
set(ASM_FILES
    src/parsing/asm.cpp
    ${BISON_ASMParser_OUTPUTS}
    ${FLEX_ASMLexer_OUTPUTS}
    )

set(SDL_FILES
    src/parsing/sdl.cpp
    ${BISON_SDLParser_OUTPUTS}
    ${FLEX_SDLLexer_OUTPUTS}
    )

add_library(qontra ${QONTRA_FILES})

if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(qontra PUBLIC -O3 -lemon -fno-strict-aliasing)
else()
    target_compile_options(qontra PUBLIC -ggdb3 -lemon -fno-strict-aliasing)
endif()

find_package(MPI REQUIRED)
find_package(LEMON REQUIRED)
find_package(Armadillo REQUIRED)

# I don't know why I have to include /usr/local/include...
target_include_directories(qontra PUBLIC
    "include" 
    ${MPI_INCLUDE_PATH} 
    "/usr/local/include"
    ${ARMADILLO_INCLUDE_DIRS})

add_subdirectory(stim)
add_subdirectory(blossom5)
add_subdirectory(cytnx)

add_executable(asm_analyzer src/main/asm_analyzer.cpp ${ASM_FILES})
add_executable(sim_test src/main/tests/sim_test.cpp ${ASM_FILES})
add_executable(cytnx_contract_test src/main/tests/cytnx/contract_test.cpp)

target_link_libraries(qontra PUBLIC libstim)
target_link_libraries(qontra PUBLIC libblossom5)
target_link_libraries(qontra PUBLIC cytnx)
target_link_libraries(qontra PUBLIC ${MPI_CXX_LIBRARIES})
target_link_libraries(qontra PUBLIC ${ARMADILLO_LIBRARIES})

target_link_libraries(asm_analyzer PRIVATE qontra)
target_link_libraries(sim_test PRIVATE qontra)
target_link_libraries(cytnx_contract_test PRIVATE qontra)
