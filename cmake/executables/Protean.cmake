# Author:   Suhas Vittal
# date:     25 December 2023

if (COMPILE_PROTEAN_LIB)
    include(${CMAKE_SOURCE_DIR}/cmake/FindCPLEX.cmake)

    set(PROTEAN_FILES
        src/protean/io.cpp
        src/protean/network/physical.cpp
        src/protean/network/raw.cpp
        src/protean/scheduler.cpp
        src/protean/experiments.cpp)
    if (PROTEAN_USE_GRAPHVIZ)
        set(PROTEAN_FILES src/protean/visualization.cpp ${PROTEAN_FILES})
    endif()
    add_library(libprotean ${PROTEAN_FILES})
    target_compile_options(libprotean PUBLIC ${COMPILE_OPTIONS})
    if (PROTEAN_PERF)
        target_compile_definitions(libprotean PUBLIC PROTEAN_PERF)
    endif()
    if (PROTEAN_USE_GRAPHVIZ)
        include(${CMAKE_SOURCE_DIR}/cmake/FindGraphviz.cmake)
        target_include_directories(libprotean PUBLIC ${GRAPHVIZ_INCLUDE_DIR})
        target_link_libraries(libprotean PUBLIC ${GRAPHVIZ_LIBRARIES})
        target_compile_definitions(libprotean PUBLIC GRAPHVIZ_ENABLED)
    endif()
    target_include_directories(libprotean PUBLIC ${CPLEX_INCLUDE_DIR})
    target_link_libraries(libprotean PUBLIC qontra 
                                            ${CPLEX_LIBRARIES})
endif()

if (COMPILE_PROTEAN_MAIN)
    add_executable(protean main/protean/main.cpp)
    target_link_libraries(protean PRIVATE libprotean)
endif()

if (COMPILE_PROTEAN_EXPERIMENTS)
    add_executable(base_memory 
        main/protean/experiments/base/memory.cpp
        src/protean/experiments.cpp)
    add_executable(planar_memory
        main/protean/experiments/base/run_planar_code.cpp)
    target_link_libraries(base_memory PRIVATE qontra)
    target_link_libraries(planar_memory PRIVATE qontra)

    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/protean)
    set_target_properties(
        base_memory planar_memory
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pr)
endif()
