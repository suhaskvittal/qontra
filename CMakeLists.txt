# Declare CMAKE version
cmake_minimum_required(VERSION 3.20.2)
#2a4e7ca86be5b659c7227ed05c2faaf0c287e807 Declare Project and Version number
project(QontraSim VERSION 0.1)
# Declare CPP Version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Declare executable.
set(QONTRA_FILES
    # Decoders
    src/decoder/mwpm.cpp
    # Parsing
    src/parsing/cmd.cpp
    # Graphs
    src/graph/decoding_graph.cpp
    src/graph/coupling_graph.cpp
    # Protean
    src/protean/compiler.cpp
    src/protean/proc3d.cpp
    src/protean/tanner_graph.cpp
    )
add_library(common INTERFACE)
add_library(qontra ${QONTRA_FILES})
add_executable(memory src/main/memory_experiment.cpp ${QONTRA_FILES})
add_executable(circuit src/main/build_circuit.cpp ${QONTRA_FILES})

if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(common INTERFACE -O3 -lemon)
else()
    target_compile_options(common INTERFACE -ggdb3 -lemon)
endif()

find_package(MPI REQUIRED)
find_package(LEMON REQUIRED)

message(${MPI_INCLUDE_PATH})

# I don't know why I have to include /usr/local/include...
target_include_directories(common INTERFACE "include" ${MPI_INCLUDE_PATH} "/usr/local/include")

add_subdirectory(stim)
add_subdirectory(blossom5)

target_link_libraries(common INTERFACE libstim)
target_link_libraries(common INTERFACE libblossom5)
target_link_libraries(common INTERFACE ${MPI_CXX_LIBRARIES})

target_link_libraries(qontra PUBLIC common)
target_link_libraries(memory PRIVATE common)
target_link_libraries(circuit PRIVATE common)
