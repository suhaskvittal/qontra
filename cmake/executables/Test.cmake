# Author:   Suhas Vittal
# date:     1 January 2024

include(CTest)

#
# PLANARITY TESTS
#

macro(make_test)
    set(SINGLE_VALUE_ARGS TARGET)
    set(MULTI_VALUE_ARGS SOURCE_FILES)
    cmake_parse_arguments(MAKE_TEST "" "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}" ${ARGN})
    message(STATUS "Compiling test ${MAKE_TEST_TARGET}...")
    
    add_executable(${MAKE_TEST_TARGET} ${MAKE_TEST_SOURCE_FILES})
    target_link_libraries(${MAKE_TEST_TARGET} PRIVATE qontra)
endmacro()

make_test(TARGET test_PLANARITY_complete SOURCE_FILES tests/graph/planar_complete.test.cpp)
make_test(TARGET test_PLANARITY_petersen SOURCE_FILES tests/graph/planar_petersen.test.cpp)

add_test(NAME PLANARITY_k4 COMMAND test_PLANARITY_complete 4)
add_test(NAME PLANARITY_k5 COMMAND test_PLANARITY_complete 5)
add_test(NAME PLANARITY_k12 COMMAND test_PLANARITY_complete 12)
add_test(NAME PLANARITY_k2_3 COMMAND test_PLANARITY_complete 2 3)
add_test(NAME PLANARITY_k3_3 COMMAND test_PLANARITY_complete 3 3)
add_test(NAME PLANARITY_k15_30 COMMAND test_PLANARITY_complete 15 30)

add_test(NAME PLANARITY_petersen COMMAND test_PLANARITY_petersen)

set_property(TEST PLANARITY_k5 PROPERTY WILL_FAIL TRUE)
set_property(TEST PLANARITY_k12 PROPERTY WILL_FAIL TRUE)
set_property(TEST PLANARITY_k3_3 PROPERTY WILL_FAIL TRUE)
set_property(TEST PLANARITY_k15_30 PROPERTY WILL_FAIL TRUE)

set_property(TEST PLANARITY_petersen PROPERTY WILL_FAIL TRUE)
