# Declare CMAKE version
cmake_minimum_required(VERSION 3.20.3)
#2a4e7ca86be5b659c7227ed05c2faaf0c287e807 Declare Project and Version number
project(Quarch VERSION 0.1)
# Declare CPP Version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Declare execution options.
option(BOOST_PATH "Path to Boost Headers" "boost_1_79_0")
option(USE_DRAMSIM3 "ON=use DRAMSim3 (default), OFF=use DRAMSim2" ON)
option(ROWHAMMER "ON=perform Rowhammer tracking" OFF)
# Declare executable.
set(QUARCH_PROGRAM_FILES
    # Main code.
    src/benchmark.cpp
    src/decoding_graph.cpp
    src/decoder.cpp
    src/mwpm_decoder.cpp
    src/defs.cpp
    src/gulliver.cpp
    # Gulliver code.
    src/gulliver/simulator.cpp
    src/gulliver/multi_qubit.cpp
    src/gulliver/cache.cpp
    src/gulliver/bdc_decoder.cpp)
add_library(quarch ${QUARCH_PROGRAM_FILES})
# Add compile options
if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(quarch PUBLIC -O3)
else()
    target_compile_options(quarch PUBLIC -ggdb3)
endif()

# Add include directory to include path.
target_include_directories(quarch PUBLIC "include")
if (BOOST_PATH)
    target_include_directories(quarch PUBLIC "${BOOST_PATH}")
else()
    # Link boost libraries to Quarch.
    find_package(boost 1.70.0 COMPONENTS graph REQUIRED)
    target_include_directories(quarch PUBLIC ${Boost_INCLUDE_DIR})
    target_link_libraries(quarch ${Boost_LIBRARIES})
endif()
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
