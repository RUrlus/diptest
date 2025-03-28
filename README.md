# diptest

[![Linux Build](https://github.com/RUrlus/diptest/actions/workflows/linux.yml/badge.svg)](https://github.com/RUrlus/diptest/actions/workflows/linux.yml)
[![Windows Build](https://github.com/RUrlus/diptest/actions/workflows/windows.yml/badge.svg)](https://github.com/RUrlus/diptest/actions/workflows/windows.yml)
[![MacOS build](https://github.com/RUrlus/diptest/actions/workflows/macos.yml/badge.svg)](https://github.com/RUrlus/diptest/actions/workflows/macos.yml)
[![PyPi](http://img.shields.io/pypi/v/diptest.svg)](https://pypi.org/project/diptest/)

A Python/C(++) implementation of Hartigan & Hartigan's dip test for unimodality.

The dip test measures multimodality in a sample by the maximum difference, over
all sample points, between the empirical distribution function, and the
unimodal distribution function that minimizes that maximum difference. Other
than unimodality, it makes no further assumptions about the form of the null
distribution.

## Usage

This library provides two functions:
* `dipstat`
* `diptest`

The first only computes Hartigan's dip statistic. `diptest` computes both the
statistic and the p-value. The p-value can be computed using interpolation of a
critical value table (default) or by bootstrapping the null hypothesis.
Note that for larger samples (N > 1e5) this is quite compute and memory intensive.

```python3
    import numpy as np
    import diptest

    # generate some bimodal random draws
    N = 1000
    hN = N // 2
    x = np.empty(N, dtype=np.float64)
    x[:hN] = np.random.normal(0.4, 1.0, hN)
    x[hN:] = np.random.normal(-0.4, 1.0, hN)

    # only the dip statistic
    dip = diptest.dipstat(x)
    
    # both the dip statistic and p-value
    dip, pval = diptest.diptest(x)
```

## Dependencies
* `numpy`
* [Optional] `OpenMP`

Parallelisation of the p-value computation using bootstrapping is offered using OpenMP.
OpenMP is disabled by default but can be enabled, see installation section
below.
 Multi-threading can be turned off by setting the number of threads equal to 1. See the docstring of `diptest` for details.


## Installation

diptest can be installed from PyPi using:

```bash
    pip install diptest
```

Wheels containing the pre-compiled extension are available for:

- Windows x84-64 - CPython 3.8 - 3.12
- Linux x84-64 - CPython 3.8 - 3.12
- MacOS x84-64 - CPython 3.8 - 3.12
- MacOS ARM-64 - CPython 3.8 - 3.12

Note that the wheels vendor/ships OpenMP with the extension to provide parallelisation out-of-the-box.
If you run into issue with multiple versions of OpenMP being loaded you have two options: build from source or install a non-bundled wheel.

### Non-bundled wheels

We provide the same wheels without OpenMP bundled here: https://github.com/RUrlus/diptest/releases 
You than install the wheel that corresponds to your Python and OS.
For example, for CPython 3.11 and MacOS ARM:

```shell
pip install diptest-0.8.0-cp311-cp311-macosx_11_0_arm64.whl
```

### Building from source

If you have a C/C++ compiler available it is advised to install without
the wheel as this enables architecture specific optimisations.

```bash
    pip install diptest --no-binary diptest
```

Compatible compilers through Pybind11:

- Clang/LLVM 3.3 or newer (for Apple Xcode's clang, this is 5.0.0 or newer)
- GCC 4.8 or newer
- Microsoft Visual Studio 2015 Update 3 or newer
- Intel classic C++ compiler 18 or newer (ICC 20.2 tested in CI)
- Cygwin/GCC (previously tested on 2.5.1)
- NVCC (CUDA 11.0 tested in CI)
- NVIDIA PGI (20.9 tested in CI)

#### 64bit index

In case you work with very large sample sizes you can use 64bit indexes.

```bash
    CMAKE_ARGS="-DDIPTEST_64BIT_index=ON" pip install diptest --no-binary diptest
```

#### Disable OpenMP

To disable OpenMP use:

```bash
    CMAKE_ARGS="-DDIPTEST_DISABLE_OPENMP=ON" pip install diptest --no-binary diptest
```

#### Debug installation

To enable a debug build use:
```bash
    CMAKE_ARGS="-DCMAKE_BUILD_TYPE=Debug" pip install diptest --no-binary diptest
```

#### Debug printing

To enable the debug print statements use:
```bash
    CMAKE_ARGS="-DDIPTEST_ENABLE_DEBUG=ON" pip install diptest --no-binary diptest
```
then call the function with debug argument set to a value greater than zero:
```python3
    diptest(x, debug=1)
```

## References

Hartigan, J. A., & Hartigan, P. M. (1985). The Dip Test of Unimodality. The
Annals of Statistics.

Hartigan, P. M. (1985). Computation of the Dip Statistic to Test for
Unimodality. Journal of the Royal Statistical Society. Series C (Applied
Statistics), 34(3), 320-325.

## Acknowledgement

`diptest` is just a Python port of [Martin Maechler's R module of the same
name](http://cran.r-project.org/web/packages/diptest/index.html).
The package wrapping the C implementation was originally written by [Alistair Muldal](https://github.com/alimuldal/diptest).
The fork is an update with a number of changes:
* Fixes a buffer overrun issue in `_dip.c` by reverting to the original C implementation
* Python bindings using Pybind11 (C++) instead of Cython
* P-value computation using bootstrapping has been moved down to C++ with optional parallelisation support through OpenMP
* Removed overhead caused by debug branching statements by placing them under a compile-time definition
* Added tests and wheel support
* C implementation of diptest was rewritten in C++ by [Prodromos Kolyvakis](https://github.com/prokolyvakis)

## License

`diptest` is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.
