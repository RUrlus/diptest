# Diptest Changelog

## v0.6.0 -- November 2023

### Enhancements

* Add full output support in diptest by @prokolyvakis in https://github.com/RUrlus/diptest/pull/29


**Full Changelog**: https://github.com/RUrlus/diptest/compare/v0.5.2...v0.6.0

## v0.5.2 -- December 2022

### Enhancements

* Added support for Python 3.11

## v0.5.1 -- June 2022

### Fix

* Fix typo in OPENMP support macro

## v0.5.0 -- June 2022

### Changes

* Disable input checks for internal calls to diptst. (Suggested by [Prodromos Kolyvakis](https://github.com/prokolyvakis))

### Enhancements

* Set `_has_open_mp_support` attribute to the extension for neater support checks

## v0.4.2 -- May 2022

### Fixes

* Fix bug in bootstrap p-value computation due to missing cast

### Changes

* Speed by moving critical values to constant class. (Special thanks to [Prodromos Kolyvakis](https://github.com/prokolyvakis))

## v0.4.1 -- May 2022

### Enhancements

* Add option to set a stream for single threaded p-value bootstrap computation

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
