# Author:   Suhas Vittal
# date:     25 December 2023

include(${CMAKE_SOURCE_DIR}/cmake/FindCPLEX.cmake)

set(PROTEAN_FILES
    src/protean/scheduler.cpp
    src/protean/utils.cpp)
add_library(protean ${PROTEAN_FILES})
target_compile_options(protean PUBLIC ${COMPILE_OPTIONS})

target_include_directories(protean PUBLIC ${CPLEX_INCLUDE_DIR})
target_link_libraries(protean PUBLIC 
                        qontra 
                        ${CPLEX_LIB_PREFIX}/libcplex.a
                        ${CPLEX_LIB_PREFIX}/libilocplex.a)

add_executable(protean_scheduler_test main/protean/tests/scheduler.cpp)
add_executable(protean_writing_test main/protean/tests/writing.cpp)

target_link_libraries(protean_scheduler_test PRIVATE protean)
target_link_libraries(protean_writing_test PRIVATE protean)
