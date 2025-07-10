/* wrapper.hpp -- header file for wrapper around diptest implementation
 * Copyright 2022 R. Urlus
 */
#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
#include <omp.h>
#endif

#include <memory>  // unique_ptr

#include <diptest/arena.hpp>
#include <diptest/common.hpp>
#include <diptest/diptest.hpp>

namespace py = pybind11;

namespace diptest {

double diptest_pval(
    const double dipstat,
    const int64_t n,
    const int64_t n_boot,
    int allow_zero,
    int debug,
    uint64_t seed,
    uint64_t stream = 0
);

double diptest_pval_mt(
    const double dipstat,
    const int64_t n,
    const int64_t n_boot,
    int allow_zero,
    int debug,
    uint64_t seed,
    size_t n_threads

);

namespace bindings {

void bind_diptest_pval(py::module& m);
void bind_diptest_pval_mt(py::module& m);

}  // namespace bindings
}  // namespace diptest
