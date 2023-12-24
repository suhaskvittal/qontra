# Copyright 2023 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.13)
project(chromobius)
include_directories(src)
set(CMAKE_CXX_STANDARD 20)

#set(SIMD_WIDTH 128)

set(CHROMOBIUS_FILES
    src/chromobius/datatypes/atomic_error.cc
    src/chromobius/datatypes/color_basis.cc
    src/chromobius/datatypes/rgb_edge.cc
    src/chromobius/decode/decoder.cc
    src/chromobius/decode/pymatcher.cc
    src/chromobius/graph/charge_graph.cc
    src/chromobius/graph/choose_rgb_reps.cc
    src/chromobius/graph/collect_atomic_errors.cc
    src/chromobius/graph/collect_composite_errors.cc
    src/chromobius/graph/collect_nodes.cc
    src/chromobius/graph/drag_graph.cc
    src/chromobius/graph/euler_tours.cc)

add_library(libchromobius ${CHROMOBIUS_FILES})


set_target_properties(libchromobius PROPERTIES PREFIX "")

target_include_directories(libchromobius PUBLIC src)
target_include_directories(libchromobius PRIVATE 
                                        ${CMAKE_SOURCE_DIR}/dependencies/stim/src
                                        ${CMAKE_SOURCE_DIR}/dependencies/pymatching/src)
target_link_libraries(libchromobius PUBLIC libstim libpym)

if(NOT(MSVC))
    target_compile_options(libchromobius PRIVATE -O3 -Wall -Wpedantic -fPIC -fno-strict-aliasing)
    target_link_options(libchromobius PRIVATE -O3)
else()
    target_compile_options(libchromobius PRIVATE)
endif()

