/* interative.hpp -- header file for iterative application of diptest
 * Copyright 2022 R. Urlus
 */
#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
#include <omp.h>
#endif

#include <algorithm>  // sort
#include <cmath>      // NAN
#include <cstdint>
#include <functional>
#include <memory>   // unique_ptr
#include <numeric>  // accumulate
#include <queue>
#include <random>     // uniform_real_distribution
#include <stdexcept>  // runtime_error
#include <vector>

#include <diptest/common.hpp>
#include <diptest/diptest.hpp>

namespace py = pybind11;

namespace diptest {

namespace details {

template <typename T>
inline void mirror_array(
    const int n,
    const int dest_size,
    const T* __restrict src,
    T* dest) {
    T* left_dest = dest;
    T* right_dest = dest + dest_size;
    // offset the dest ptr st we are in the middle
    T* left_center = left_dest + (n - 1);
    T* right_center = right_dest + (n - 1);
    T* right_end = right_dest + dest_size;
    // the center is always zero
    *left_center = 0.0;
    *right_center = 0.0;

    const T min = *std::min_element(src, src + n);
    const T max = *std::max_element(src, src + n);
    for (int li = 1, ri = 0; li < n; ++li, ++ri) {
        T left_val = src[li] - min;
        // mirror left
        left_center[-li] = -left_val;
        left_center[li] = left_val;

        // mirror right
        T right_val = src[ri] - max;
        right_dest[ri] = right_val;
        right_end[-li] = std::abs(right_val);
    }
}

}  // namespace details

template <typename T>
inline py::array_t<T> mirror_array(const py::array_t<T>& src) {
    const int src_size = src.size();
    const int dest_size = (2 * src_size - 1);

    auto dest = py::array_t<T, py::array::f_style>({dest_size, 2});
    details::mirror_array(src_size, dest_size, src.data(), dest.mutable_data());
    return dest;
}

py::array_t<double> iterdip(
    const py::array_t<double>& x,
    const double crit_val,
    const int n_boot,
    const int max_iter,
    const int allow_zero,
    const int debug,
    uint64_t seed,
    uint64_t stream);

namespace bindings {

void bind_mirror_array(py::module& m);
void bind_iterdip(py::module& m);

}  // namespace bindings
}  // namespace diptest
