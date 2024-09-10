# author:   Suhas Vittal
# date:     6 January 2024
#
# Basic utilities library.

cmake_minimum_required(VERSION 3.20.2)
project(Vtils VERSION 0.1)

if (CMAKE_BUILD_TYPE MATCHES Release)
    set(COMPILE_OPTIONS -O3 -fno-strict-aliasing)
else()
    set(COMPILE_OPTIONS -ggdb3 -Wall)
endif()
                        
set(CMAKE_CXX_STANDARD 20)

set(VTILS_SOURCE_FILES 
    src/vtils/cmd_parse.cpp
    src/vtils/filesystem.cpp
    src/vtils/ini_parse.cpp
    src/vtils/mat2.cpp
    src/vtils/timer.cpp
    src/vtils/linprog/base.cpp)

add_library(vtils ${VTILS_SOURCE_FILES})
target_compile_options(vtils PRIVATE ${COMPILE_OPTIONS})
target_include_directories(vtils PUBLIC "include")

