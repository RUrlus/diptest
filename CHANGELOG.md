# Diptest Changelog

## 0.10.0 -- July 2025

### Features

* FEAT: Add dip index/location to full output by @JohannesBuchner  in https://github.com/RUrlus/diptest/pull/78

### CICD

* CICD: Bump pypa/cibuildwheel from 2.23.2 to 2.23.3 by @dependabot in https://github.com/RUrlus/diptest/pull/73
* CICD: Bump pypa/cibuildwheel from 2.23.3 to 3.0.0 by @dependabot in https://github.com/RUrlus/diptest/pull/74
* CICD: Bump pypa/cibuildwheel from 3.0.0 to 3.0.1 by @dependabot in https://github.com/RUrlus/diptest/pull/75

**Full Changelog**: https://github.com/RUrlus/diptest/compare/v0.9.0...v0.10.0


## 0.9.0 -- April 2025

### Features

* FEAT: Offer compile-time option to use 64bit indexes

## 0.8.3 -- March 2025

### Fixes

* FIX: Resolve negative dip statistic values due to overflow by @RUrlus in https://github.com/RUrlus/diptest/pull/69
* FIX: Correct flag check for contiguous memory layout by @RUrlus in https://github.com/RUrlus/diptest/pull/69

### CICD

* CICD: Bump Clang and GCC versions by @RUrlus in https://github.com/RUrlus/diptest/pull/67
* CICD: Bump pypa/cibuildwheel from 2.20.0 to 2.21.0 by @dependabot in https://github.com/RUrlus/diptest/pull/59
* CICD: Bump pypa/cibuildwheel from 2.21.0 to 2.21.1 by @dependabot in https://github.com/RUrlus/diptest/pull/60
* CICD: Bump pypa/cibuildwheel from 2.21.1 to 2.21.2 by @dependabot in https://github.com/RUrlus/diptest/pull/61
* CICD: Bump pypa/cibuildwheel from 2.21.3 to 2.22.0 by @dependabot in https://github.com/RUrlus/diptest/pull/64
* CICD: Bump pypa/cibuildwheel from 2.22.0 to 2.23.1 by @dependabot in https://github.com/RUrlus/diptest/pull/66

## v0.8.2 -- November 2024

* BLD: use relative path for libomp by @RUrlus in https://github.com/RUrlus/diptest/pull/63

## v0.8.1 -- August 2024

* PKG: Suppress error when two version of OpenMP are loaded

## v0.8.0 -- April 2024

* CHG: [C++] Use the cheap-multiplier variant of the PCG64DXSM by @RUrlus in https://github.com/RUrlus/diptest/pull/46
* BLD: [C++] Replace outdated architecture flag script and build-out CMake files by @RUrlus in https://github.com/RUrlus/diptest/pull/47
* BLD: Vendor OpenMP binary in wheel by @RUrlus in https://github.com/RUrlus/diptest/pull/39

### CICD

* Bump actions/setup-python from 3 to 5 by @dependabot in https://github.com/RUrlus/diptest/pull/41
* Bump pypa/cibuildwheel from 2.16.2 to 2.16.5 by @dependabot in https://github.com/RUrlus/diptest/pull/40
* Bump pypa/gh-action-pypi-publish from 1.8.11 to 1.8.14 by @dependabot in https://github.com/RUrlus/diptest/pull/43
* Bump pypa/cibuildwheel from 2.16.5 to 2.17.0 by @dependabot in https://github.com/RUrlus/diptest/pull/44

## v0.7.0 -- December 2023

* ENH: Add support for Python 3.12 by @rurlus
* MAINT: Switch to scikit-build-core as build system by @rurlus
* FIX: Correct off-by-one error in the indexes of the full results by @rurlus

## v0.6.1 -- November 2023

### Fix

* FIX: Correct the out of index error in interpolation by @prokolyvakis in https://github.com/RUrlus/diptest/pull/32

## v0.6.0 -- November 2023

### Enhancements

* Add full output support in diptest by @prokolyvakis in https://github.com/RUrlus/diptest/pull/29

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
