# Grid

## for the SU(2) Nf=1 and Nf=2 DiRAC project

**Data parallel C++ mathematical object library.**

License: GPL v2.

Last update June 2017.

For the latest version of the underlying Grid code, please see [the main repository][Grid]

## Building

This was built on Tursa using CUDA 11.4.1, Open MPI 4.1.1, and UCX 1.12.0.

```
module load cuda/11.4.1 openmpi/4.1.1-cuda11.4.1  ucx/1.12.0-cuda11.4.1
```

To rebuild:

```
git clone https://github.com/edbennett/Grid_epcc Grid
cd Grid
mkdir build
cd build
../configure \
    --enable-comms=mpi \
    --enable-simd=GPU \
    --enable-shm=nvlink \
    --enable-gen-simd-width=64 \
    --enable-accelerator=cuda \
    --enable-Nc=2 \
    --enable-accelerator-cshift \
    --disable-unified \
    --disable-gparity\
    CXX=nvcc \
    LDFLAGS="-cudart shared" \
    CXXFLAGS="-ccbin mpicxx -gencode arch=compute_80,code=sm_80 -std=c++14 -cudart shared"
make -j24
cd tests/hmc
make Test_hmc_EOWilsonFermionGauge Test_rhmc_Wilson1p1
```

## Running

### Executables
* For Nf=1: use Test_rhmc_Wilson1p1
* For Nf=2: use Test_hmc_EOWilsonFermionGauge

### Optimal parallelisations by lattice volume

* L=24: 1 GPU. Four jobs can be executed in parallel on a 4-GPU node.
  * If insufficient runs are available (e.g. before runs thermalise), 2 GPUs can be used
    without loss of speed.
* L=32: 8 GPUs, two nodes.
* L=48: 32 GPUs, eight nodes.

Example batch scripts (and corresponding wrapper scripts) for single and multi-node jobs are in the `runscripts` directory.

### Initialisation

Unlike HiRep, Grid gives very poor acceptance after a hot start. To overcome this, set `--Thermalizations 20` on the first execution.
This should be removed before subsequent continuation runs when restarting from a checkpoint.

### Multiple streams

At any point after thermalisation a stream can be split to give two parallel streams, allowing more rapid results given sufficient resources.

To do this, set some random seeds, e.g. using `--serialseed "11 12 13 14 15" --parallelseed "16 17 18 19 20"`, and use `--StartingType CheckpointStartWithReseed`
to tell Grid to reseed the generators with these. (The default seeds are `--serialseed "1 2 3 4 5" --parallelseed "6 7 8 9 10"`.)


[Grid]: https://github.com/paboyle/Grid