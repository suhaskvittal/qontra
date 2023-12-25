# Copyright 2021 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.13)
project(stim)
include_directories(src)
set(CMAKE_CXX_STANDARD 20)

# Convert desired SIMD_WIDTH into machine architecture flags.
if (ARCH_x86_64)
    if(SIMD_WIDTH EQUAL 256)
        set(MACHINE_FLAG "-mavx2" "-msse2")
    elseif(SIMD_WIDTH EQUAL 128)
        set(MACHINE_FLAG "-mno-avx2" "-msse2")
    elseif(SIMD_WIDTH EQUAL 64)
        set(MACHINE_FLAG "-march=native -mno-avx2" "-mno-sse2")
    endif()
else()
    set(MACHINE_FLAG "")
endif()

# make changes to file_lists trigger a reconfigure
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS file_lists/source_files_no_main)

file(STRINGS file_lists/source_files_no_main SOURCE_FILES_NO_MAIN)

add_library(libstim ${SOURCE_FILES_NO_MAIN})
set_target_properties(libstim PROPERTIES PREFIX "")
target_include_directories(libstim PUBLIC src)

target_compile_definitions(libstim PUBLIC SIMD_WIDTH=${SIMD_WIDTH})

if(NOT(MSVC))
    target_compile_options(libstim PRIVATE -O3 -Wall -Wpedantic -fPIC -fno-strict-aliasing ${MACHINE_FLAG})
    target_link_options(libstim PRIVATE -O3)
else()
    target_compile_options(libstim PRIVATE ${MACHINE_FLAG})
endif()
