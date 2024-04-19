# Author:   Suhas Vittal
# date:     25 December 2023

macro(make_main)
    set(SINGLE_VALUE_ARGS TARGET)
    set(MULTI_VALUE_ARGS SOURCE_FILES)
    cmake_parse_arguments(MAKE_MAIN "" "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}" ${ARGN})

    add_executable(${MAKE_MAIN_TARGET} ${MAKE_MAIN_SOURCE_FILES})
    target_link_libraries(${MAKE_MAIN_TARGET} PRIVATE qontra)
endmacro()

make_main(TARGET converter SOURCE_FILES main/converter.cpp)
make_main(TARGET generate_syndromes SOURCE_FILES main/generate_syndromes.cpp)
make_main(TARGET memory SOURCE_FILES main/memory.cpp)
make_main(TARGET qontrasim SOURCE_FILES main/qontrasim.cpp)
make_main(TARGET print_error_stats SOURCE_FILES main/print_error_stats.cpp)

make_main(TARGET explain_degenerate_errors SOURCE_FILES main/explain_degenerate_errors.cpp)
make_main(TARGET compute_code_distance SOURCE_FILES main/compute_code_distance.cpp)
make_main(TARGET compute_circuit_distance SOURCE_FILES main/compute_circuit_distance.cpp)
make_main(TARGET decoder_debugger SOURCE_FILES main/decoder_debugger.cpp)
