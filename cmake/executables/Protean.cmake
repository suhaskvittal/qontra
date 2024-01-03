# Author:   Suhas Vittal
# date:     25 December 2023

include(${CMAKE_SOURCE_DIR}/cmake/FindCPLEX.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FindGraphviz.cmake)

set(PROTEAN_FILES
    src/protean/network/io.cpp
    src/protean/network/physical.cpp
    src/protean/network/raw.cpp
    src/protean/visualization.cpp
    src/protean/visualization_attr.cpp
    src/protean/visualization_lp.cpp)
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
