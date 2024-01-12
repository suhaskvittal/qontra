# Author:   Suhas Vittal
# date:     25 December 2023

# If the dependencies are packaged with Qontra, then add them here:

add_subdirectory(dependencies/stim)
add_subdirectory(dependencies/blossom5)
add_subdirectory(dependencies/qes)
add_subdirectory(dependencies/vtils)

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
target_link_libraries(qontra PUBLIC qes)
target_link_libraries(qontra PUBLIC vtils)

target_link_libraries(qontra PUBLIC ${MPI_CXX_LIBRARIES})

if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(qontra PUBLIC "stdc++fs")
endif()

if (COMPILE_NEURAL_DECODER)
    target_include_directories(qontra PUBLIC ${ARMADILLO_INCLUDE_DIRS} ${MLPACK_INCLUDE_DIRS})
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
