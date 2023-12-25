# Author:   Suhas Vittal
# date:     25 December 2023

# Setup compile options for Qontra and executables.
if (NOT SIMD_WIDTH)
    set(SIMD_WIDTH 64)
endif()

if (NOT ARCH_x86_64)
    if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
        set(ARCH_x86_64 ON)
    endif()
endif()

if (NOT COMPILE_OPTIONS) 
    set(RELEASE_COMPILE_OPTIONS -Ofast -fno-strict-aliasing -Wno-psabi)
    set(DEBUG_COMPILE_OPTIONS -ggdb3 -fno-strict-aliasing -Wall)

    if (CMAKE_BUILD_TYPE MATCHES "Release")
        set(COMPILE_OPTIONS ${RELEASE_COMPILE_OPTIONS})
    else()
        set(COMPILE_OPTIONS ${DEBUG_COMPILE_OPTIONS})
    endif()
endif()

message(STATUS "(Build type = ${CMAKE_BUILD_TYPE}, compile options are: ${COMPILE_OPTIONS}")

# Set CXX standard. Default is C++17, but certain dependencies (if enabled) require C++20.
set(CMAKE_CXX_STANDARD 17)

if (COMPILE_PYMATCHING)
    set(CMAKE_CXX_STANDARD 20)
endif()

if (COMPILE_CHROMOBIUS)
    set(CMAKE_CXX_STANDARD 20)
    set(COMPILE_PYMATCHING On)  # PyMatching is a dependency.
endif()

set(CMAKE_CXX_STANDARD_REQUIRED True)
