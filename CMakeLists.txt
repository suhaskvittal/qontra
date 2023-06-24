# Declare CMAKE version
cmake_minimum_required(VERSION 3.20.2)
#2a4e7ca86be5b659c7227ed05c2faaf0c287e807 Declare Project and Version number
project(QontraSim VERSION 0.1)
# Declare CPP Version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Declare executable.
set(QONTRA_FILES
    # Top-Level
    src/experiments.cpp
    # Decoders
    src/decoder/mwpm.cpp
    # Parsing
    src/parsing/cmd.cpp
    # Graphs
    src/graph/decoding_graph.cpp src/graph/coupling_graph.cpp
    src/graph/tanner_graph.cpp
    # Protean
    src/protean/compiler.cpp
    src/protean/proc3d.cpp
    # Simulation
    src/sim/clifford_sim.cpp
    src/sim/control_sim.cpp
    )
add_library(qontra ${QONTRA_FILES})
add_executable(memory src/main/memory_experiment.cpp)
add_executable(circuit src/main/build_circuit.cpp)

### TESTING ###
add_executable(cliff src/main/test/cliffsimtest.cpp)
add_executable(ctrl src/main/test/ctrlsimtest.cpp)

### Protean ###
add_executable(protean_main src/main/protean/main.cpp)

if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(qontra PUBLIC -O3 -lemon -fno-strict-aliasing)
else()
    target_compile_options(qontra PUBLIC -ggdb3 -lemon -fno-strict-aliasing)
endif()

find_package(MPI REQUIRED)
find_package(LEMON REQUIRED)

message(${MPI_INCLUDE_PATH})

# I don't know why I have to include /usr/local/include...
target_include_directories(qontra PUBLIC "include" ${MPI_INCLUDE_PATH} "/usr/local/include")

add_subdirectory(stim)
add_subdirectory(blossom5)

target_link_libraries(qontra PUBLIC libstim)
target_link_libraries(qontra PUBLIC libblossom5)
target_link_libraries(qontra PUBLIC ${MPI_CXX_LIBRARIES})

target_link_libraries(memory PRIVATE qontra)
target_link_libraries(circuit PRIVATE qontra)
target_link_libraries(cliff PRIVATE qontra)
target_link_libraries(ctrl PRIVATE qontra)

target_link_libraries(protean_main PRIVATE qontra)
