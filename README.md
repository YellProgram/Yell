
Yell: a program for diffuse scattering interpretation
=======

Yell is a program for analyzing diffuse scattering from single crystals using Three Dimensional Difference Pair Distribution Function  (3D-∆PDF) method.

### Executable binary files
You can download the latest version of Yell for Mac and Windows [here](https://github.com/YellProgram/Yell/releases/latest).

### Documentation
The pdf version of documentation available [here](https://github.com/YellProgram/Yell/releases/download/v1.1/Yell.reference.pdf).

### Examples
[The hypothetical iron-void example](https://arkadiysimonov.github.io/yellfiles/iron-void_example.zip).
The resulting PDF files can be visualized using the [PDFViewer](https://github.com/aglie/DensityViewer/releases) or using data analysis languages like [python](https://arkadiysimonov.github.io/yellfiles/Visualize_yell_python.zip).

---

## Building from source

The build downloads several external dependencies (Eigen, HDF5, bitshuffle, glog, gflags,
Ceres Solver) via CMake `ExternalProject`. Because these downloads happen at build time,
running a parallel build (`-j16`) on a clean directory can race — `yell-lib` sources may
start compiling before HDF5/bitshuffle headers are written to disk. Use the sequences
below to avoid this.

### Prerequisites

- CMake ≥ 3.10
- C++20-capable compiler (GCC 10+, Clang 13+, MSVC 2022+)
- Boost 1.88 unpacked at `../libs/boost_1_88_0` relative to the source root
- On macOS: Xcode command-line tools (provides Accelerate framework for LAPACK)
- On Linux: `lapack`, `blas`, `gfortran`

> **Forcing static HDF5 download**: By default, if a system HDF5 is found it will be
> used. To force downloading and building HDF5 statically (recommended for reproducible
> builds), add `-DCMAKE_TOTAL_STATIC=ON` to the configure step.

### Debug build

```bash
cd /path/to/Yell

# 1. Configure (add -DCMAKE_TOTAL_STATIC=ON to force static HDF5 download)
cmake -B cmake-build-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_TOTAL_STATIC=ON

# 2. Download and build HDF5 (must finish before bitshuffle can configure)
cmake --build cmake-build-debug --target download_hdf5 -j1

# 3. Build bitshuffle (needs HDF5 headers from step 2)
cmake --build cmake-build-debug --target build_bitshuffle -j1

# 4. Download Eigen (needed for yell-lib compilation)
cmake --build cmake-build-debug --target eigen -j16

# 5. Build glog, gflags, then Ceres (Ceres depends on both)
cmake --build cmake-build-debug --target glog-lib --target extern_gflags -j16
cmake --build cmake-build-debug --target ceres-solver -j16

# 6. Build the main binary
cmake --build cmake-build-debug --target yell -j16
```

### Release build

```bash
cd /path/to/Yell

# 1. Configure (add -DCMAKE_TOTAL_STATIC=ON to force static HDF5 download)
cmake -B cmake-build-release -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOTAL_STATIC=ON

# 2. Download and build HDF5
cmake --build cmake-build-release --target download_hdf5 -j1

# 3. Build bitshuffle
cmake --build cmake-build-release --target build_bitshuffle -j1

# 4. Download Eigen
cmake --build cmake-build-release --target eigen -j16

# 5. Build glog, gflags, then Ceres
cmake --build cmake-build-release --target glog-lib --target extern_gflags -j16
cmake --build cmake-build-release --target ceres-solver -j16

# 6. Build the main binary
cmake --build cmake-build-release --target yell -j16
```

### Subsequent builds (after first download)

Once all external projects have been downloaded and built, CMake skips them on
subsequent runs. You can rebuild just the main binary quickly:

```bash
cmake --build cmake-build-debug --target yell -j16   # debug
cmake --build cmake-build-release --target yell -j16  # release
```

### If something goes wrong

To force a re-download of HDF5/bitshuffle without wiping the whole build directory:

```bash
rm -rf cmake-build-debug/hdf5 cmake-build-debug/download_bitshuffle-prefix cmake-build-debug/build_bitshuffle-prefix
# then redo steps 2–6 above
```


