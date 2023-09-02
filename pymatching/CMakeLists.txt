cmake_minimum_required(VERSION 3.4...3.18)
project(PyMatching2)

include_directories(src)
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(SIMD_WIDTH 128)

set(ARCH_OPT "-O2")
# https://stackoverflow.com/a/8591946 (RTC1 incompatible with O2 flag)
STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

set(SOURCE_FILES_NO_MAIN
        src/pymatching/sparse_blossom/driver/io.cc
        src/pymatching/sparse_blossom/driver/mwpm_decoding.cc
        src/pymatching/sparse_blossom/flooder/graph.cc
        src/pymatching/sparse_blossom/flooder/detector_node.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.cc
        src/pymatching/sparse_blossom/flooder/graph_fill_region.cc
        src/pymatching/sparse_blossom/flooder/match.cc
        src/pymatching/sparse_blossom/flooder/graph_flooder.cc
        src/pymatching/sparse_blossom/matcher/alternating_tree.cc
        src/pymatching/sparse_blossom/matcher/mwpm.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/region_edge.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.cc
        src/pymatching/sparse_blossom/tracker/flood_check_event.cc
        src/pymatching/sparse_blossom/search/search_graph.cc
        src/pymatching/sparse_blossom/search/search_detector_node.cc
        src/pymatching/sparse_blossom/search/search_flooder.cc
        src/pymatching/sparse_blossom/driver/user_graph.cc
        src/pymatching/rand/rand_gen.cc
        )

add_library(libpym ${SOURCE_FILES_NO_MAIN})

target_compile_options(libpym PRIVATE ${ARCH_OPT})
if (NOT (MSVC))
    target_link_options(libpym PRIVATE -pthread ${ARCH_OPT})
endif ()
target_link_libraries(libpym PUBLIC libstim)
target_include_directories(libpym PUBLIC src)
