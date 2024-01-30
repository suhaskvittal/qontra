# Author:   Suhas Vittal
# date:     25 August 2022

cmake_minimum_required(VERSION 3.20.2)
project(QontraSim VERSION 0.2)

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

if (EXISTS ${CMAKE_SOURCE_DIR}/cmake/UserConfig.cmake)
    message(STATUS "Loading variables from user's config.")
    include(${CMAKE_SOURCE_DIR}/cmake/UserConfig.cmake)
endif()

set(QONTRA_ISA_FILE "${CMAKE_SOURCE_DIR}/isa.txt")

include(${CMAKE_SOURCE_DIR}/cmake/ConfigureOptions.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/SetupQontra.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/LinkDependencies.cmake)

include(${CMAKE_SOURCE_DIR}/cmake/executables/Main.cmake)


if (COMPILE_PROTEAN)
    include(${CMAKE_SOURCE_DIR}/cmake/executables/Protean.cmake)
endif()

if (COMPILE_PROBIUS)
    include(${CMAKE_SOURCE_DIR}/cmake/executables/Probius.cmake)
endif()

if (COMPILE_TESTS)
    include(${CMAKE_SOURCE_DIR}/cmake/executables/Test.cmake)
endif()

if (COMPILE_SSE)
    include(${CMAKE_SOURCE_DIR}/cmake/executables/SSE.cmake)
endif()
