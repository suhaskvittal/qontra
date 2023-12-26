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

# Make changes to file_lists trigger a reconfigure.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS file_lists/source_files_no_main)
file(STRINGS file_lists/source_files_no_main SOURCE_FILES_NO_MAIN)

add_library(libchromobius ${SOURCE_FILES_NO_MAIN})
set_target_properties(libchromobius PROPERTIES PREFIX "")
target_include_directories(libchromobius PUBLIC src 
                            ${CMAKE_SOURCE_DIR}/dependencies/stim/src
                            ${CMAKE_SOURCE_DIR}/dependencies/pymatching/src)
target_link_libraries(libchromobius PRIVATE libstim libpym)
if(NOT(MSVC))
    target_compile_options(libchromobius PRIVATE -O3 -Wall -Wpedantic -fPIC -fno-strict-aliasing)
    target_link_options(libchromobius PRIVATE -O3)
else()
    target_compile_options(libchromobius PRIVATE)
endif()
