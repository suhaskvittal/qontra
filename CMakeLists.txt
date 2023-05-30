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
    src/graph/decoding_graph.cpp)
add_library(qontra ${QONTRA_FILES})
add_executable(memory src/main/memory_experiment.cpp ${QONTRA_FILES})
add_executable(circuit src/main/build_circuit.cpp ${QONTRA_FILES})
# Add compile options
if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(qontra PUBLIC -O3)
    target_compile_options(memory PUBLIC -O3)
else()
    target_compile_options(qontra PUBLIC -ggdb3)
    target_compile_options(memory PUBLIC -ggdb3)
endif()
target_compile_options(circuit PUBLIC -O3)

find_package(MPI REQUIRED)

# Add include directory to include path.
target_include_directories(qontra PUBLIC "include" ${MPI_INCLUDE_PATH})
target_include_directories(memory PUBLIC "include" ${MPI_INCLUDE_PATH})
target_include_directories(circuit PUBLIC "include" ${MPI_INCLUDE_PATH})

add_subdirectory(stim)
target_link_libraries(qontra libstim)
target_link_libraries(memory libstim)
target_link_libraries(circuit libstim)

add_subdirectory(blossom5)
target_link_libraries(qontra libblossom5)
target_link_libraries(memory libblossom5)
target_link_libraries(circuit libblossom5)

#add_subdirectory(dramsim3)
#target_link_libraries(quarch dramsim3)

target_link_libraries(qontra ${MPI_CXX_LIBRARIES})
target_link_libraries(memory ${MPI_CXX_LIBRARIES})
