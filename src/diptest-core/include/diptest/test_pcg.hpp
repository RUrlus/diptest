/* test_pcg.hpp -- header file for tests of PCG RNG as used in wrapper.cpp
 * Copyright 2022 R. Urlus */

#ifndef INCLUDE_DIPTEST_TEST_PCG_HPP_
#define INCLUDE_DIPTEST_TEST_PCG_HPP_
#if defined(DIPTEST_BUILD_CPP_TESTS)

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
#include <omp.h>
#endif

#include <cstddef>
#include <random>  // uniform_real_distribution

#include <diptest/common.hpp>

namespace py = pybind11;

namespace diptest {
namespace tests {

py::array_t<double> pcg_seed_test(const size_t size, const uint64_t seed);
py::array_t<double> pcg_set_stream_test(
    const size_t size,
    const uint64_t seed,
    const py::array_t<uint64_t>& streams);
#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
py::array_t<double> pcg_mt_stream_test(
    const size_t row_size,
    const uint64_t seed,
    const int n_threads);
#endif  // DIPTEST_HAS_OPENMP_SUPPORT

}  // namespace tests

namespace bindings {
void bind_pcg_seed_test(py::module& m);
void bind_pcg_set_stream_at_init(py::module& m);
void bind_pcg_set_stream_test(py::module& m);
#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
void bind_pcg_mt_stream_test(py::module& m);
#endif  // DIPTEST_HAS_OPENMP_SUPPORT
}  // namespace bindings
}  // namespace diptest

#endif  // DIPTEST_BUILD_CPP_TESTS
#endif  // INCLUDE_DIPTEST_TEST_PCG_HPP_
