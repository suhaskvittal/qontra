# Top-Level Files and Directories

This README details the contents of the top-level files and their directories.
For more information about each directory, go to the directory's README.

## File structure:

```
qontra
|-- decoder.h and decoder 
|-- experiments.h and experiments 
|-- ext 
|-- graph.h and graph
|-- sim 
|-- tables.h
```

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

## Graphs

Unfortunately, to the best of our knowledge, there is no "good" C++ graph
library. Initially, QontraSim used Boost, but the extreme compile times for
Boost Graph coupled with odd segfaults motivated a custom graph library.
QontraSim's graph library boasts the following features:
1. Low feature bloat: QontraSim's graph library isn't intended to be the Swiss army
   knife that many other graph libraries aspire to be. For instance, other
   libraries allow fine-grained optimizations such as changing backing data
   structures. QontraSim, for the sake of *simplicity*, offers no such
   optimizations. This simplicity avoids the problem of obtuse syntax as seen
   with Boost Graph.
2. Memory-safe: all vertices/edges are implemented via smart pointers.
   It is very hard to cause a segfault.
3. Extensible: the basic graph class is templated so that common use cases can
   be fulfilled by simply defining vertices and edges. Furthermore, if the user
   requires more control, then they can simply extend the base class.

### Vertices and Edges

`graph.h` implements the base class `Graph` along with other objects. We first
discuss the concept of vertices and edges:
```c++
namespace base {

struct vertex_t {
    uint64_t id;
};

struct edge_t {
    sptr<void> src;
    sptr<void> dst;
    bool is_undirected=true;

    template <class V=void> sptr<V> get_source(void);
    template <class V=void> sptr<V> get_target(void);
};

}
```
As claimed before, these classes are rather barebones. `base::vertex_t` has 
a unique `id` that is used to identify a vertex within a Graph. `base::edge_t`
has two endpoints and a boolean indicating if the edge is undirected. If a user
wants to add additional information to the vertices or edges, they must
implement a new class. It is *recommended* that new vertex and edge classes
extend `base::vertex_t` and `base::edge_t`, respectively.

The oddest implementation detail here is that `base::edge_t` implements its
endpoints as `void` pointers. This is done to avoid implicit dependencies
between `base::vertex_t` and `base::edge_t`. Instead, the user can use whatever
vertex they want with `base::edge_t` and vice versa. Consequently, if the user
needs an endpoint, they can retrieve it safely with `get_source` or
`get_target`.

Vertices and edges in general also have the following utility functions for
printing information.
```c++
template <class V> std::string                  print_v(sptr<V>);
template <> std::string                         print_v(sptr<void>);
template <class V=void, class E> std::string    print_e(sptr<E>);
```
If the user defines a new vertex or edge class and would like to modify the
output information, these functions can be specialized for those classes.

### Graphs

The base class `Graph` is a barebones graph class that offers many common
functions. The following is a list of these functions:
```c++
virtual void    change_id(sptr<V>, uint64_t to);
virtual void    manual_update_id(sptr<V>, uint64_t old_id, uint64_t new_id);

virtual bool    contains(uint64_t id);
virtual bool    contains(sptr<V>);
virtual bool    contains(sptr<V>, sptr<V>);
virtual bool    contains(sptr<E>);

virtual sptr<V> make_vertex(void);
virtual sptr<E> make_edge(sptr<V>, sptr<V>, bool is_undirected=true);

virtual bool    add_vertex(sptr<V>);
virtual bool    add_edge(sptr<E>);

virtual sptr<V> make_and_add_vertex(uint64_t id);
virtual sptr<E> make_and_add_edge(sptr<V>, sptr<V>, bool is_undirected=true);

virtual sptr<V> get_vertex(uint64_t);
virtual sptr<E> get_edge(sptr<V>, sptr<V>);
virtual sptr<E> get_edge(uint64_t, uint64_t);

virtual void    delete_vertex(sptr<V>);
virtual void    delete_edge(sptr<E>);

std::vector<sptr<V>>    get_vertices(void);
std::vector<sptr<E>>    get_edges(void);

size_t  n(void);    // Returns number of vertices.
size_t  m(void);    // Returns number of edges.

std::vector<sptr<V>>    get_neighbors(sptr<V>);
std::vector<sptr<V>>    get_incoming(sptr<V>);
std::vector<sptr<V>>    get_outgoing(sptr<V>);

std::vector<sptr<V>>    get_common_neighbors(sptr<V>, sptr<V>);

size_t  get_degree(sptr<V>);
size_t  get_indegree(sptr<V>);
size_t  get_outdegree(sptr<V>);
size_t  get_inoutdegree(sptr<V>);

fp_t    get_mean_degree(void);
size_t  get_max_degree(void);
```
Our implementation of `Graph` uses an adjacency matrix for O(1) edge lookups, an
adjacency list for O(1) neighborhood retrieval, and a reverse adjacency list to
track incoming edges in digraphs. Most functions are O(1), except for the
`delete_vertex` and `delete_edge` functions, which are O(n) and O(m),
respectively.

`Graph` also has the capability to maintain a sense of state. This is rather
useful in cases where a property of a graph is repeatedly accessed but is
expensive to recompute (i.e. planarity). Thus, for subclasses of `Graph`, the
protected function
```c++
virtual bool update_state(void);
```
can be overridden to track custom properties. `update_state` returns `true` if
the `Graph` was updated, and `false` otherwise. Thus, when overriding this
function, it is recommended that the function body follows the format:
```c++
bool update_state() {
    if (!<PARENT>::update_state()) return false;
    // Else update state:
    ...
    return true;
}
```
where `<PARENT>` is the immediate superclass. The implementation of
`update_state` in `Graph` obviously does not follow the above format, but
instead only updates the state of the graph if the protected variable
`graph_has_changed` is `true`. This variable is set on modifications, such as
vertex/edge additions or deletions. If subclasses also need to indicate when the
graph has changed, it is recommended to either update `graph_has_changed` or use
their own variable. Furthermore, the public function
```c++
void force_update_state(void);
```
is also provided so that users can force an update, regardless of the value of
`graph_has_changed`.

## Simulations

## Error Modeling
