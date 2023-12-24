cmake_minimum_required(VERSION 3.4...3.18)
project(PyMatching2)

include_directories(src)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCE_FILES_NO_MAIN
        src/pymatching/sparse_blossom/driver/namespaced_main.cc
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
        src/pymatching/sparse_blossom/diagram/animation_main.cc
        src/pymatching/sparse_blossom/diagram/mwpm_diagram.cc
        src/pymatching/sparse_blossom/search/search_graph.cc
        src/pymatching/sparse_blossom/search/search_detector_node.cc
        src/pymatching/sparse_blossom/search/search_flooder.cc
        src/pymatching/sparse_blossom/driver/user_graph.cc
        src/pymatching/rand/rand_gen.cc
        )

add_library(libpym ${SOURCE_FILES_NO_MAIN})
set_target_properties(libpym PROPERTIES PREFIX "")
target_include_directories(libpym PUBLIC src)
target_compile_options(libpym PRIVATE -O3 -fPIC)
if(NOT(MSVC))
    target_link_options(libpym PRIVATE -pthread -O3)
endif()
target_link_libraries(libpym libstim)
