# Author:   Suhas Vittal
# date:     1 January 2024

include(CTest)

macro(make_test)
    set(SINGLE_VALUE_ARGS TARGET)
    set(MULTI_VALUE_ARGS SOURCE_FILES)
    cmake_parse_arguments(MAKE_TEST "" "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}" ${ARGN})
    message(STATUS "Compiling test ${MAKE_TEST_TARGET}...")
    
    add_executable(${MAKE_TEST_TARGET} ${MAKE_TEST_SOURCE_FILES})
    target_link_libraries(${MAKE_TEST_TARGET} PRIVATE qontra)
endmacro()

make_test(TARGET test_PLANARITY SOURCE_FILES tests/graph/planar_complete.test.cpp)

add_test(NAME test_PLANARITY_K4 COMMAND test_PLANARITY 4)
add_test(NAME test_PLANARITY_K5 COMMAND test_PLANARITY 5)
add_test(NAME test_PLANARITY_K6 COMMAND test_PLANARITY 6)

add_test(NAME test_PLANARITY_K2_3 COMMAND test_PLANARITY 2 3)
add_test(NAME test_PLANARITY_K3_3 COMMAND test_PLANARITY 3 3)

set_property(TEST test_PLANARITY_K5 PROPERTY WILL_FAIL TRUE)
set_property(TEST test_PLANARITY_K6 PROPERTY WILL_FAIL TRUE)
set_property(TEST test_PLANARITY_K3_3 PROPERTY WILL_FAIL TRUE)
