# QontraSim v0.2

QontraSim is a general-purpose simulator for Quantum Error Correction, with a
focus on the control system. QontraSim is partially built upon Stim, and boasts
the following features:
1. A library of decoders, including matching decoders and neural network
   decoders.
2. Circuit simulators, such as re-implementation of Stim's `FrameSimulator` and
   a Clifford simulator.
3. A graph library.

## Dependencies

Most of QontraSim's dependencies are packaged with the source code and are
modified versions of the original code (i.e. Stim and PyMatching). Other
dependencies are:
1. CMake v3.20.2 and C++20 are required to compile QontraSim. QontraSim has been
   tested with GCC *only*.
2. OpenMPI, which can be installed via package managers (i.e. homebrew or apt).
   OpenMPI is used by QontraSim to parallelize simulations.
3. [MLPACK](https://www.mlpack.org/getstarted.html) for neural network
   decoders. MLPACK should be installed by the user, and the install directory
   should be provided via the Cmake argument `MLPACK_INCLUDE_DIRS` (see
   [here](#building-qontrasim)).
4. [QES (Quantum Experiment Specification)](https://github.com/suhaskvittal/qes)
   is used generally by QontraSim as a way of specifying experiments. QES is
   packaged with QontraSim as a submodule.
5. [Vtils](https://github.com/suhaskvittal/vtils) is a set of utility functions
   that are used by QontraSim. Vtils is also packaged with QontraSim as a
   submodule.

## Building QontraSim

QontraSim uses CMake as the build system. For most installations, the following
commands will build a Release version of QontraSim.
```
$ mkdir build && cd build
$ cmake .. -DCMAKE_BUILD_TYPE=Release
$ make -j4
```
The user can also specify different options in `cmake/UserConfig.cmake`, which
is loaded at the beginning of the build. Here, the user can set the following
options:
* `COMPILE_NEURAL_DECODER`: if this option is set, then the user must provide
  `MLPACK_INCLUDE_DIRS` and have MLPACK and its dependencies installed on their
  system.
* `COMPILE_PYMATCHING` and `COMPILE_CHROMOBIUS`: if these options are set, then
  wrappers for [PyMatching](https://github.com/oscarhiggott/PyMatching) and
  [Chromobius](https://github.com/quantumlib/chromobius) are made available.
* `L1D_CACHE_LINE_SIZE`: this option is used to optimize the batch size for
  cache locality. The size must be given in **bytes**, and is available via
  `lscpu` on Linux or `sysctl -a | grep cachelinesize` on MacOS.

### Other Notes

Because QontraSim mostly uses integer and bitwise operations, the Release build
is compiled with `-Ofast`; when `COMPILE_NEURAL_DECODER` is set, this is set to
`-O3`. If the user wishes to use `-O3` instead of `-Ofast`, then
`RELEASE_COMPILE_OPTIONS` should be modified in `cmake/ConfigureOptions.cmake`.

## Compiled Executables

QontraSim comes with *four* executables for general use.
1. `converter`, which converts an input file into an output file. This is useful
   for debugging `.qes` files by converting them into `.stim` files.
2. `generate_syndromes`, which spits out syndromes generated by running the
   instructions in a `.qes` file. This is also rather useful for debugging and
   verifying if syndrome extraction circuits are functionally correct.
3. `memory`, which runs a memory experiment given a `.qes` file using a matching
   decoder.
4. `qontrasim`, which runs a full system simulation given a `.qes` file. The
   outputs are (1) a folder of syndromes in the `.dets` format (aka "traces"),
   (2) a probability histogram printed to stdout, (3) a `.stim` file that can be
   used by a decoder to decode the given syndromes, and (4) additional data in
   an output file.

