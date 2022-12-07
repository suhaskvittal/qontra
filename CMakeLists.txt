# Declare CMAKE version
cmake_minimum_required(VERSION 3.20.3)
#2a4e7ca86be5b659c7227ed05c2faaf0c287e807 Declare Project and Version number
project(Quarch VERSION 0.1)
# Declare CPP Version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Declare executable.
set(QUARCH_PROGRAM_FILES
    # Main code.
    src/benchmark.cpp
    src/decoding_graph.cpp
    src/decoder.cpp
    src/mwpm_decoder.cpp
    src/defs.cpp
    # HYPERION
    src/hyperion.cpp
    src/hyperion/simulator.cpp
    src/hyperion/multi_qubit.cpp
    src/hyperion/trial_decoder.cpp
    # TMR
    src/tmr_decoder.cpp
    # FLEECE
    src/fleece.cpp
    src/fleece/lattice_graph.cpp)
add_library(quarch ${QUARCH_PROGRAM_FILES})
# Add compile options
if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(quarch PUBLIC -O3)
else()
    target_compile_options(quarch PUBLIC -ggdb3)
endif()

# Add include directory to include path.
target_include_directories(quarch PUBLIC "include")
# Link OpenMP to Quarch
if (ENABLE_OPENMP)
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        target_link_libraries(quarch OpenMP::OpenMP_CXX)
    endif()
endif()
# Add Stim to Quarch.
add_subdirectory(stim)
target_link_libraries(quarch libstim)
# Add Blossom V to Quarch
add_subdirectory(blossom5)
target_link_libraries(quarch libblossom5)
# Add DRAMSim to Quarch
add_subdirectory(dramsim3)
target_link_libraries(quarch dramsim3)
