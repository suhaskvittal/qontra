# Author:   Suhas Vittal
# date:     25 December 2023

if (COMPILE_PROTEAN_MAIN)
    include(${CMAKE_SOURCE_DIR}/cmake/FindCPLEX.cmake)
    include(${CMAKE_SOURCE_DIR}/cmake/FindGraphviz.cmake)

    set(PROTEAN_FILES
        src/qontra/protean/io.cpp
        src/qontra/protean/network/physical.cpp
        src/qontra/protean/network/raw.cpp
        src/qontra/protean/scheduler.cpp
        src/qontra/protean/visualization.cpp
        src/qontra/protean/visualization_attr.cpp
        src/qontra/protean/visualization_lp.cpp
        src/qontra/protean/yield_sim.cpp
        src/qontra/protean/experiments.cpp)
    add_library(libprotean ${PROTEAN_FILES})
    target_compile_options(libprotean PUBLIC ${COMPILE_OPTIONS})
#   target_compile_definitions(libprotean PUBLIC PROTEAN_PERF)
    target_include_directories(libprotean PUBLIC ${CPLEX_INCLUDE_DIR}
                                                    ${GRAPHVIZ_INCLUDE_DIR})
    target_link_libraries(libprotean PUBLIC qontra 
                                            ${CPLEX_LIBRARIES}
                                            ${GRAPHVIZ_LIBRARIES})

    add_executable(protean main/protean/main.cpp)
    target_link_libraries(protean PRIVATE libprotean)

    target_compile_definitions(protean PRIVATE GRAPHVIZ_ENABLED)
endif()

if (COMPILE_PROTEAN_EXPERIMENTS)
    add_executable(pr_base_memory 
        main/protean/experiments/base/memory.cpp
        src/qontra/protean/experiments.cpp)
    add_executable(pr_planar_memory
        main/protean/experiments/base/run_planar_code.cpp)
    target_link_libraries(pr_base_memory PRIVATE qontra)
    target_link_libraries(pr_planar_memory PRIVATE qontra)

#   target_compile_definitions(pr_base_memory PUBLIC PROTEAN_PERF)
endif()
