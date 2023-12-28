/* test_pcg.cpp -- implementation of testss for PCG RNG as used in wrapper.cpp
 * Copyright 2022 R. Urlus */
#if defined(DIPTEST_BUILD_CPP_TESTS)
#include <diptest/test_pcg.hpp>

namespace py = pybind11;

namespace diptest {
namespace tests {

py::array_t<double> pcg_seed_test(const size_t size, const uint64_t seed) {
    // allocate the return array
    py::array_t<double> arr(size);
    // declare and initialise the RNG
    details::pcg64_dxsm rng;
    if (seed == 0) {
        details::pcg_seed_seq seed_source;
        rng.seed(seed_source);
    } else {
        rng.seed(seed);
    }
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double* ptr = arr.mutable_data();
    for (size_t i = 0; i < size; i++) {
        ptr[i] = dist(rng);
    }
    return arr;
}

py::array_t<double> pcg_set_stream_test(
    const size_t size,
    const uint64_t seed,
    const py::array_t<uint64_t>& streams) {
    const size_t n_streams = streams.size();
    // allocate the return array
    py::array_t<double> arr(size * n_streams);
    // declare and initialise the RNG
    details::pcg64_dxsm rng;
    if (seed == 0) {
        details::pcg_seed_seq seed_source;
        rng.seed(seed_source);
    } else {
        rng.seed(seed);
    }
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double* ptr = arr.mutable_data();
    for (size_t j = 0; j < n_streams; j++) {
        // offset ptr
        // set stream
        rng.set_stream(streams.at(j));
        for (size_t i = 0; i < size; i++) {
            ptr[i] = dist(rng);
        }
        ptr += size;
    }
    return arr.reshape({n_streams, size});
}

#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
py::array_t<double> pcg_mt_stream_test(
    const size_t row_size,
    const uint64_t seed,
    const int n_threads) {
    details::pcg64_dxsm global_rng;
    if (seed == 0) {
        details::pcg_seed_seq seed_source;
        global_rng.seed(seed_source);
    } else {
        global_rng.seed(seed);
    }

    const size_t total_size = row_size * n_threads;
    py::array_t<double> arr(total_size);
    double* ptr = arr.mutable_data();

#pragma omp parallel num_threads(n_threads) shared(ptr, global_rng)
    {
        details::pcg64_dxsm rng = global_rng;
        rng.set_stream(omp_get_thread_num() + 1);
        std::uniform_real_distribution<double> dist(0.0, 1.0);

#pragma omp for
        for (size_t i = 0; i < total_size; i++) {
            ptr[i] = dist(rng);
        }
    }  // pragma parallel
    return arr.reshape({static_cast<size_t>(n_threads), row_size});
}
#endif  // DIPTEST_HAS_OPENMP_SUPPORT

}  // namespace tests

namespace bindings {

void bind_pcg_seed_test(py::module& m) {
    m.def(
        "pcg_seed_test",
        &diptest::tests::pcg_seed_test,
        py::arg("size"),
        py::arg("seed"));
}

void bind_pcg_set_stream_test(py::module& m) {
    m.def(
        "pcg_set_stream_test",
        &diptest::tests::pcg_set_stream_test,
        py::arg("size"),
        py::arg("seed"),
        py::arg("streams"));
}

#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
void bind_pcg_mt_stream_test(py::module& m) {
    m.def(
        "pcg_mt_stream_test",
        &diptest::tests::pcg_mt_stream_test,
        py::arg("row_size"),
        py::arg("seed"),
        py::arg("n_threads"));
}
#endif

}  // namespace bindings
}  // namespace diptest

#endif  // DIPTEST_BULD_CPP_TESTS
