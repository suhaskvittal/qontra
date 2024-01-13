# Author:   Suhas Vittal
# date:     25 December 2023

include(${CMAKE_SOURCE_DIR}/cmake/FindCPLEX.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FindGraphviz.cmake)

set(PROTEAN_FILES
    src/qontra/protean/io.cpp
    src/qontra/protean/network/physical.cpp
    src/qontra/protean/network/raw.cpp
    src/qontra/protean/scheduler.cpp
    src/qontra/protean/visualization.cpp
    src/qontra/protean/visualization_attr.cpp
    src/qontra/protean/visualization_lp.cpp)
add_library(libprotean ${PROTEAN_FILES})
target_compile_options(libprotean PUBLIC ${COMPILE_OPTIONS})

target_include_directories(libprotean PUBLIC ${CPLEX_INCLUDE_DIR}
                                                ${GRAPHVIZ_INCLUDE_DIR})
target_link_libraries(libprotean PUBLIC qontra 
                                        ${CPLEX_LIBRARIES}
                                        ${GRAPHVIZ_LIBRARIES})

add_executable(protean main/protean/main.cpp)
target_link_libraries(protean PRIVATE libprotean)

target_compile_definitions(protean PRIVATE GRAPHVIZ_ENABLED)

if (COMPILE_PROTEAN_EXPERIMENTS)
    add_executable(pr_nn_train main/protean/experiments/nn_train.cpp)
    target_link_libraries(pr_nn_train PRIVATE qontra)

    add_executable(pr_nn_memory main/protean/experiments/nn_memory.cpp)
    target_link_libraries(pr_nn_memory PRIVATE qontra)
endif()
