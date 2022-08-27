# Declare CMAKE version
cmake_minimum_required(VERSION 3.23.2)
# Declare compiler
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)
# Declare Project and Version number
project(Quarch VERSION 0.1)
# Declare CPP Version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Declare execution options.
option(ENABLE_OPENMP "use openmp" OFF)
option(MWPM_QUANTIZE_LEVEL "mwpm quantization level (0=None, 1=16b, 2=8b)" 0)
option(USE_DRAMSIM3 "ON=use DRAMSim3 (default), OFF=use DRAMSim2" ON)
# Declare executable.
set(QUARCH_PROGRAM_FILES
    # Main code.
    src/main.cpp
    src/stim_test.cpp
    src/benchmark.cpp
    src/decoding_graph.cpp
    src/decoder.cpp
    src/mwpm_decoder.cpp
    src/defs.cpp
    src/lilliput.cpp
    src/gulliver.cpp
    src/clique.cpp
    # Gulliver code.
    src/gulliver/experiments.cpp)
add_executable(quarch ${QUARCH_PROGRAM_FILES})
# Add compile options
if (ENABLE_OPENMP) 
    target_compile_options(quarch PRIVATE -O3 -Wall -fopenmp)
    target_compile_definitions(quarch PRIVATE USE_OMP)
elseif(CMAKE_BUILD_TYPE EQUAL Release)
    target_compile_options(quarch PRIVATE -O3 -Wall)
endif()
target_compile_definitions(quarch PRIVATE 
                            QFP_SIZE=${MWPM_QUANTIZE_LEVEL}
                            HOME_DIRECTORY="${CMAKE_SOURCE_DIR}")
# Add include directory to include path.
target_include_directories(quarch PRIVATE "include")
# Link boost libraries to Quarch.
find_package(boost 1.79.0 COMPONENTS graph REQUIRED)
target_include_directories(quarch PUBLIC ${Boost_INCLUDE_DIR})
target_link_libraries(quarch ${Boost_LIBRARIES})
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
if (USE_DRAMSIM3)
    add_subdirectory(dramsim3)
    target_link_libraries(quarch dramsim3)
else()
    add_subdirectory(dramsim2)
    target_link_libraries(quarch libdramsim2)
endif()
