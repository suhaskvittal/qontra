# Top-Level Files and Directories

This README details the contents of the top-level files and their directories.
For more information about each directory, go to the directory's README.

## File structure:

qontra

|-- [decoder.h and decoder](#decoders) 

|-- [experiments.h and experiments](#experiments) 

|-- [ext](#qontra-extensions) 

|-- [graph.h and graph](#graphs) 
j
|-- [sim](#simulation) 

|-- [tables.h](#error-modeling) 

## Decoders

Decoders are constructed as subclasses of the virtual `Decoder` class. The
`Decoder` constructor is of the form:
```c++
Decoder(DetailedStimCircuit, graph::DecodingGraph::Mode);
```
`DetailedStimCircuit` is an extension of Stim's `Circuit` class (for more
information, look [here](#qontra-extensions)), but supports casting from
`stim::Circuit` (so `stim::Circuit` can be passed directly for the first
argument). `graph::DecodingGraph::Mode` specifies how the decoder will build the
decoding graph for the underlying circuit. *By default*, this is set to
`graph::DecodingGraph::Mode::NORMAL`, which constructs the distance matrix for
the decoding graph. Note that for large circuits, this is rather memory
intensive, and not all decoders may need the distance matrix. To avoid
constructing the distance matrix, use `graph::DecodingGraph::Mode::LOW_MEMORY`.
If the decoder does not need the decoding graph at all (i.e. neural nets), use
`graph::DecodingGraph::Mode::DO_NOT_BUILD`.

To implement a decoding algorithm, a decoder must implement:
```c++
Decoder::result_t Decoder::decode_error(stim::simd_bits_range_ref<SIMD_WIDTH>);
```
Note that `SIMD_WIDTH` is a constant depending on the CMake configuration.
Nevertheless, the `Decoder::result_t` contains three entries, not all of which
need to be set:
1. `fp_t exec_time`, or the amount of time taken (in nanoseconds) to decode a
   syndrome. Note that `Decoder` offers a protected member `vtils::Timer timer`
   to subclasses for timing purposes. To start the timer, use
   `timer.clk_start()`. To end timing, use `timer.clk_end()`, which returns a
   `uint64_t` as the time elapsed in nanoseconds.
2. `stim::simd_bits<SIMD_WIDTH> corr`, the correction for the given syndrome.
   The width of this bitvector is equal to the number of observables for the
   given circuit.
3. `std::vector<assign_t> error_assignments`, or the matching-based assignments
   made by the decoder. As this is inherently specific to matching-based
   decoders, this is not guaranteed to be set.

Finally, to retrieve the nonzero detectors from a syndrome, `decoder.h` offers
the `get_nonzero_detectors_(T, uint64_t)` and `Decoder::get_nonzero_detectors(T)`,
the latter of which is just `get_nonzero_detectors_(T, circuit.count_detectors())`.
The first argument, which is a template class, should either be `stim::simd_bits` or
`stim::simd_bits_range_ref`. For `get_nonzero_detectors_`, the second argument
is the number of detectors in syndrome.

## Experiments

Experiments, such as *syndrome generation* or *memory experiments*, can be run
via functions in `experiments.h` and other files in `experiments`.
`experiments.h` contains the basic functions for running experiments, and
the functions for other experiments generally stem from these functions.

`experiments.h` has five global parameters that can be set by the user:
1. `G_USE_MPI`: if set to false, then MPI is not used. Note that the **user**
   must call `MPI_Init` and `MPI_Finalize`. Default is `true`.
2. `G_SHOTS_PER_BATCH`: number of shots performed at once. The default is
   `100'000`, but this can be configured a bit more optimally by calling
   `configure_optimal_batch_size`, which sets `G_SHOTS_PER_BATCH` to the number
   of bits in a cacheline (i.e. 512 bits on most processors, but some
   processors, like Apple's Mx series may have 1 Kbit cachelines).
3. `G_BASE_SEED`: important when using MPI. To avoid having all processes
   generate the same syndromes, the seed of the *k*-th process is set to
   `G_BASE_SEED + k`. Default is `0`.
4. `G_FILTER_OUT_SYNDROMES`: ignore syndromes with sufficiently low Hamming
   weight. This is useful to get quick results for memory experiments when low
   Hamming weight syndromes are unlikely to cause errors. Default is `true`.
5. `G_FILTERING_HAMMING_WEIGHT`: if `G_FILTER_OUT_SYNDROMES` is set, then all
   syndromes with Hamming weights less than or equal to
   `G_FILTERING_HAMMING_WEIGHT` are ignored. By default, this is 0.

As mentioned above in point (2), `configure_optimal_batch_size` can be called to
fit the batch size to the number of bits in the L1d cacheline. On Linux, this
can be found by calling `sysconf(_SC_LEVEL1_DCACHE_LINESIZE)`, but for systems
where this does not work, the user should set `L1D_CACHE_LINE_SIZE` when using
CMake. Note that `L1D_CACHE_LINE_SIZE` is the size in *bytes*, not bits; by
default, this is set to 64 bytes.

`experiments.h` also has three functions for generating and performing I/O with
syndromes:
```c++
// This is the input to all callbacks shown here.
struct shot_payload_t {
    stim::simd_bits_range_ref<SIMD_WIDTH> syndrome;
    stim::simd_bits_range_ref<SIMD_WIDTH> observables;
};

void        generate_syndromes(const DetailedStimCircuit&, uint64_t shots, CALLBACK);
uint64_t    read_syndrome_trace(std::string input_folder, const DetailedStimCircuit&, CALLBACK);
void        build_syndrome_trace(std::string output_folder, const DetailedStimCircuit&, uint64_t shots);
```
As with `Decoder`, `stim::Circuit` can be used in place of `DetailedStimCircuit`
without issues. `generate_syndromes` is straightforward and generates as many
syndromes as shots specified. The user can further specify a callback
(`CALLBACK`) that takes in a `shot_payload_t` that can be used to interact with
the generated syndromes. `read_syndrome_trace` operates similarly, but with the
syndromes in a folder containing `.dets` files (see Stim file formats). Finally,
`build_syndrome_trace` generates syndromes and writes them to a folder
containing `.dets` files. Note that `build_syndrome_trace` makes one `.dets`
file *per* batch so that reading the syndromes can be quickly parallelized when
calling `read_syndrome_trace`.

## Qontra Extensions

QontraSim provides extensions to the Stim's circuit model. The first is moreso
for custom simulations supported by QontraSim via the [QES (Quantum Experiment
Specification)](https://github.com/suhaskvittal/qes) language, which is built
for QontraSim. QES is intended to be modular and useful for writing arbitrary
experiments without having to make significant modifications to simulation
infrastructure nor the specific error rates. The second is the
`DetailedStimCircuit` class, which extends `stim::Circuit` by adding additional
information, such as indicating which qubits are flag qubits. We believe this is
more readable and extensible than Stim's current strategy of handling extra
information, which encoding this additional information into the coordinates of
the detectors.

# 
