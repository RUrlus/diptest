# Diptest Changelog

## v0.4.1 -- May 2022

### Enhancements

* Add option to set a stream for single threaded p-value bootstrap computation

### Changes 

* diptest.c was rewritten in C++ (Special thanks to [Prodromos Kolyvakis](https://github.com/prokolyvakis))
* Incorporated OptimizeForArchitecture from VC for better architecture specific
  compile flags

## v0.4.0 -- May 2022

### Changes 

* diptest.c was rewritten in C++ (Special thanks to [Prodromos Kolyvakis](https://github.com/prokolyvakis))
* Incorporated OptimizeForArchitecture from VC for better architecture specific
  compile flags

## v0.3.0 -- April 2022

### Changes

* Switch to PCG64-DXSM RNG from Mersenne twister

## v0.2.3 -- April 2022

Patch release

### Changes

* Fix conversion to double in accumulate

## v0.2.2 -- March 2022

Patch release

### Changes

* Fix for incorrect number of default threads in bootstrap p-value computation
* Minimal scikit-build version is 0.14.1

#### Internal

* Reduce memory footprint single-threaded bootstrap computation p-value

## v0.2.1 -- March 2022

Patch release

### Changes

* Enforce C99 standard in CMake

## 0.2.0 -- March 2022

Initial release of the fork of https://github.com/alimuldal/diptest

### Changes

* Fixes a buffer overrun issue in `_dip.c` by reverting to the original C implementation
* Python bindings using Pybind11 (C++) instead of Cython

### Enhancements

* P-value computation using bootstrapping has been moved down to C++ with optional parallelisation support through OpenMP
* Removed overhead caused by debug branching statements by placing them under a compile-time definition
* Added tests and wheel support
