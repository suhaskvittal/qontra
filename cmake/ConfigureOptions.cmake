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
    set(RELEASE_COMPILE_OPTIONS -fno-strict-aliasing -Wno-psabi -Wno-write-strings)
    set(DEBUG_COMPILE_OPTIONS -ggdb3 -fno-strict-aliasing -Wall -Wno-write-strings -Wno-bool-operation)

    if (CMAKE_BUILD_TYPE MATCHES "Release")
        if (COMPILE_NEURAL_DECODER)
            set(CXX_OPTIMIZATION_LEVEL -O3)
        else()
            set(CXX_OPTIMIZATION_LEVEL -Ofast)
        endif()
        set(COMPILE_OPTIONS ${CXX_OPTIMIZATION_LEVEL} ${RELEASE_COMPILE_OPTIONS})
    else()
        set(COMPILE_OPTIONS ${DEBUG_COMPILE_OPTIONS})
    endif()
endif()

message(STATUS "(Build type = ${CMAKE_BUILD_TYPE}, compile options are: ${COMPILE_OPTIONS})")

set(CMAKE_CXX_STANDARD 20)

if (COMPILE_CHROMOBIUS)
    set(COMPILE_PYMATCHING On)  # PyMatching is a dependency.
endif()

set(CMAKE_CXX_STANDARD_REQUIRED True)
