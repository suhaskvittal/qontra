# Author:   Suhas Vittal
# date:     25 December 2023

# If the dependencies are packaged with Qontra, then add them here:

add_subdirectory(dependencies/stim)
add_subdirectory(dependencies/blossom5)

if (COMPILE_PYMATCHING)
    message(STATUS "Will compile and link PyMatching.")
    add_subdirectory(dependencies/pymatching)
endif()

if (COMPILE_CHROMOBIUS)
    message(STATUS "Will compile and link Chromobius.")
    add_subdirectory(dependencies/chromobius)
endif()

# Add include directories and link Qontra to dependencies.

target_include_directories(qontra PUBLIC "include" ${MPI_INCLUDE_PATH})

target_link_libraries(qontra PUBLIC libstim)
target_link_libraries(qontra PUBLIC libblossom5)
target_link_libraries(qontra PUBLIC ${MPI_CXX_LIBRARIES})

if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(qontra PUBLIC "stdc++fs")
endif()

if (COMPILE_NEURAL_DECODER)
    target_include_directories(qontra PUBLIC ${ARMADILLO_INCLUDE_DIRS})
    target_include_directories(qontra PUBLIC ${ARMADILLO_INCLUDE_DIRS})
    target_link_libraries(qontra PUBLIC ${ARMADILLO_LIBRARIES})

    if (OpenMP_CXX_FOUND)
        message(STATUS "Found OpenMP for the NN Decoder.")
        target_link_libraries(qontra PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()

if (COMPILE_PYMATCHING)
    target_link_libraries(qontra PUBLIC libpym)
endif()

if (COMPILE_CHROMOBIUS)
    target_link_libraries(qontra PUBLIC libchromobius)
endif()

if (LINK_LEMON)
    if (DEFINED LEMON_INSTALL_DIR)
        set(LEMON_LIBRARY_DIR "${LEMON_INSTALL_DIR}/lemon/lib")
        set(LEMON_INCLUDE_DIR "${LEMON_INSTALL_DIR}/lemon/include")
    endif()

    if (NOT DEFINED LEMON_LIBRARY_DIR)
        message(FATAL_ERROR "COMPILE_LEMON is set, but LEMON_LIBRARY_DIR is not defined.")
    elseif(NOT DEFINED LEMON_INCLUDE_DIR)
        message(FATAL_ERROR "COMPILE_LEMON is set, but LEMON_INCLUDE_DIR is not defined.")
    else()
        # Everything checks out. Link lemon.
        message(STATUS "Will link Lemon, LEMON_LIBRARY_DIR=${LEMON_LIBRARY_DIR}, LEMON_INCLUDE_DIR=${LEMON_INCLUDE_DIR}")
        add_library(lemon STATIC IMPORTED)
        set_target_properties(lemon PROPERTIES IMPORTED_LOCATION "${LEMON_LIBRARY_DIR}/libemon.a")
        target_link_libraries(qontra PUBLIC lemon)
        target_include_directories(qontra PUBLIC ${LEMON_INCLUDE_DIR})
    endif()
endif()
