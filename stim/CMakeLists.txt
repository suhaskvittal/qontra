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
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY out)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY out)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY out)

# Convert desired SIMD_WIDTH into machine architecture flags.
if(NOT(SIMD_WIDTH))
    set(MACHINE_FLAG "")
#   set(MACHINE_FLAG "-march=native")
elseif(SIMD_WIDTH EQUAL 256)
    set(MACHINE_FLAG "-mavx2" "-msse2")
elseif(SIMD_WIDTH EQUAL 128)
    set(MACHINE_FLAG "-mno-avx2" "-msse2")
elseif(SIMD_WIDTH EQUAL 64)
    set(MACHINE_FLAG "-mno-avx2" "-mno-sse2")
endif()

set(SOURCE_FILES_NO_MAIN
        src/stim/arg_parse.cc
        src/stim/circuit/circuit.cc
        src/stim/circuit/gate_data.cc
        src/stim/circuit/gate_data_annotations.cc
        src/stim/circuit/gate_data_blocks.cc
        src/stim/circuit/gate_data_collapsing.cc
        src/stim/circuit/gate_data_controlled.cc
        src/stim/circuit/gate_data_hada.cc
        src/stim/circuit/gate_data_noisy.cc
        src/stim/circuit/gate_data_pauli.cc
        src/stim/circuit/gate_data_period_3.cc
        src/stim/circuit/gate_data_period_4.cc
        src/stim/circuit/gate_data_pp.cc
        src/stim/circuit/gate_data_swaps.cc
        src/stim/circuit/gate_target.cc
        src/stim/dem/detector_error_model.cc
        src/stim/help.cc
        src/stim/gen/circuit_gen_main.cc
        src/stim/gen/circuit_gen_params.cc
        src/stim/gen/gen_color_code.cc
        src/stim/gen/gen_rep_code.cc
        src/stim/gen/gen_surface_code.cc
        src/stim/io/measure_record_batch.cc
        src/stim/io/measure_record_batch_writer.cc
        src/stim/io/measure_record.cc
        src/stim/io/measure_record_reader.cc
        src/stim/io/measure_record_writer.cc
        src/stim/io/raii_file.cc
        src/stim/io/stim_data_formats.cc
        src/stim/io/sparse_shot.cc
        src/stim/main_namespaced.cc
        src/stim/mem/bit_ref.cc
        src/stim/mem/simd_bit_table.cc
        src/stim/mem/simd_bits.cc
        src/stim/mem/simd_bits_range_ref.cc
        src/stim/mem/simd_compat.cc
        src/stim/mem/simd_util.cc
        src/stim/mem/sparse_xor_vec.cc
        src/stim/probability_util.cc
        src/stim/search/graphlike/algo.cc
        src/stim/search/graphlike/edge.cc
        src/stim/search/graphlike/graph.cc
        src/stim/search/graphlike/node.cc
        src/stim/search/graphlike/search_state.cc
        src/stim/search/hyper/algo.cc
        src/stim/search/hyper/edge.cc
        src/stim/search/hyper/graph.cc
        src/stim/search/hyper/node.cc
        src/stim/search/hyper/search_state.cc
        src/stim/simulators/dem_sampler.cc
        src/stim/simulators/detection_simulator.cc
        src/stim/simulators/error_analyzer.cc
        src/stim/simulators/error_matcher.cc
        src/stim/simulators/frame_simulator.cc
        src/stim/simulators/matched_error.cc
        src/stim/simulators/measurements_to_detection_events.cc
        src/stim/simulators/tableau_simulator.cc
        src/stim/simulators/vector_simulator.cc
        src/stim/stabilizers/conversions.cc
        src/stim/stabilizers/pauli_string.cc
        src/stim/stabilizers/pauli_string_ref.cc
        src/stim/stabilizers/tableau.cc
        src/stim/stabilizers/tableau_iter.cc
        src/stim/stabilizers/tableau_specialized_prepend.cc
        src/stim/stabilizers/tableau_transposed_raii.cc
        src/stim.cc
        )

add_library(libstim ${SOURCE_FILES_NO_MAIN})
set_target_properties(libstim PROPERTIES PREFIX "")
target_include_directories(libstim PUBLIC src)
if(NOT(MSVC))
    target_compile_options(libstim PRIVATE -O3 -Wall -Wpedantic -fno-strict-aliasing)
    target_link_options(libstim PUBLIC -pthread -O3)
endif()
target_compile_options(libstim PUBLIC ${MACHINE_FLAG})
