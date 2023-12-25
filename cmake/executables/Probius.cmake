# Author:   Suhas Vittal
# date:     25 December 2023

# There are executable-specific options because PyMatching needs C++20, but LEMON does not work with
# C++20 for whatever reason.
if (COMPILE_PROBIUS_TRACER_ONLY)
    set(COMPILE_PROBIUS_TRACER ON)
    set(COMPILE_PROBIUS_MEMORY OFF)
elseif (COMPILE_PROBIUS_MEMORY_ONLY)
    set(COMPILE_PROBIUS_TRACER OFF)
    set(COMPILE_PROBIUS_MEMORY ON)
else()
    set(COMPILE_PROBIUS_TRACER ON)
    set(COMPILE_PROBIUS_MEMORY ON)
endif()

if (COMPILE_PROBIUS_TRACER)
    add_executable(probius_tracer main/probius/tracer.cpp)
    target_link_libraries(probius_tracer PRIVATE qontra)
    target_compile_options(probius_tracer PRIVATE ${COMPILE_OPTIONS})
endif()

if (COMPILE_PROBIUS_MEMORY)
    add_executable(probius_memory main/probius/memory.cpp)
    target_link_libraries(probius_memory PRIVATE ${COMPILE_OPTIONS})
endif()
