/* wrapper.cpp -- implementation of wrapper around diptst from diptest.c
 * Copyright 2022 R. Urlus
 */
#include <diptest/wrapper.hpp>

namespace py = pybind11;

namespace diptest {

namespace details {

inline double diptest(const double* x_ptr, int N, int allow_zero, int debug) {
    std::array<int, 4> lo_hi = {0, 0, 0, 0};
    std::unique_ptr<int[]> gcm(new int[N]);
    std::unique_ptr<int[]> lcm(new int[N]);
    std::unique_ptr<int[]> mn(new int[N]);
    std::unique_ptr<int[]> mj(new int[N]);

    double dip = diptst<true>(
        x_ptr,
        N,
        lo_hi.data(),
        gcm.get(),
        lcm.get(),
        mn.get(),
        mj.get(),
        allow_zero,
        debug);
    return dip;
}  // diptest

}  // namespace details

double diptest(const py::array_t<double>& x, int allow_zero, int debug) {
    return details::diptest(x.data(), x.size(), allow_zero, debug);
}  // diptest

py::dict diptest_full(const py::array_t<double>& x, int allow_zero, int debug) {
    const double* x_ptr = x.data();
    int N = x.size();
    std::array<int, 4> lo_hi = {0, 0, 0, 0};

    auto gcm = py::array_t<int>(N);
    auto lcm = py::array_t<int>(N);
    int* gcm_ptr = gcm.mutable_data();
    int* lcm_ptr = lcm.mutable_data();

    std::unique_ptr<int[]> mn(new int[N]);
    std::unique_ptr<int[]> mj(new int[N]);
    int* mn_ptr = mn.get();
    int* mj_ptr = mj.get();

    double dip = diptst<true>(
        x_ptr,
        N,
        lo_hi.data(),
        gcm_ptr,
        lcm_ptr,
        mn_ptr,
        mj_ptr,
        allow_zero,
        debug);

    using namespace pybind11::literals;  // to bring in the `_a` literal NOLINT
    return py::dict(
        "dip"_a = dip,
        "lo"_a = lo_hi[0],
        "hi"_a = lo_hi[1],
        "xl"_a = x.at(lo_hi[0]),
        "xu"_a = x.at(lo_hi[1]),
        "_gcm"_a = gcm,
        "_lcm"_a = lcm,
        "_lh_2"_a = lo_hi[2],
        "_lh_3"_a = lo_hi[3]);
}  // diptest_full

double diptest_pval(
    const double dipstat,
    const int64_t n,
    const int64_t n_boot,
    int allow_zero,
    int debug,
    uint64_t seed,
    uint64_t stream) {
    pcg64_dxsm rng;
    if (seed == 0) {
        pcg_seed_seq seed_source;
        rng.seed(seed_source);
    } else if (stream != 0) {
        rng.seed(seed, stream);
    } else {
        rng.seed(seed);
    }
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double dip;
    std::array<int, 4> lo_hi = {0, 0, 0, 0};
    std::unique_ptr<int[]> gcm(new int[n]);
    std::unique_ptr<int[]> lcm(new int[n]);
    std::unique_ptr<int[]> mn(new int[n]);
    std::unique_ptr<int[]> mj(new int[n]);
    std::unique_ptr<int[]> dips(new int[n_boot]);
    std::unique_ptr<double[]> sample(new double[n]);

    double* r_sample = sample.get();
    double* sample_end = r_sample + n;

    for (int64_t i = 0; i < n_boot; i++) {
        for (int64_t j = 0; j < n; j++) {
            r_sample[j] = dist(rng);
        }
        std::sort(r_sample, sample_end);
        dip = diptst<false>(
            r_sample,
            n,
            lo_hi.data(),
            gcm.get(),
            lcm.get(),
            mn.get(),
            mj.get(),
            allow_zero,
            debug);
        dips[i] = dipstat <= dip;
    }
    int64_t accu = 0;
    double p_val = static_cast<double>(
                       std::accumulate(dips.get(), dips.get() + n_boot, accu))
                   / n_boot;
    return p_val;
}  // diptest_pval

#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
double diptest_pval_mt(
    const double dipstat,
    const int64_t n,
    const int64_t n_boot,
    int allow_zero,
    int debug,
    uint64_t seed,
    size_t n_threads) {
    std::unique_ptr<bool[]> dips(new bool[n_boot]);
    pcg64_dxsm global_rng;
    if (seed == 0) {
        pcg_seed_seq seed_source;
        global_rng.seed(seed_source);
    } else {
        global_rng.seed(seed);
    }

#pragma omp parallel num_threads(n_threads) shared(dips, global_rng)
    {
        std::unique_ptr<int[]> lo_hi(new int[n]);
        std::memset(lo_hi.get(), 0, 4);
        std::unique_ptr<int[]> gcm(new int[n]);
        std::unique_ptr<int[]> lcm(new int[n]);
        std::unique_ptr<int[]> mn(new int[n]);
        std::unique_ptr<int[]> mj(new int[n]);
        std::unique_ptr<double[]> sample(new double[n]);

        double* p_sample = sample.get();
        double* p_sample_end = p_sample + n;

        // PCG family has different streams which are, in theory, independent of
        // each other. Hence, we can use the same seed and a different stream to
        // draw independent samples from each thread without having to allocate
        // the whole block
        pcg64_dxsm rng = global_rng;
        rng.set_stream(omp_get_thread_num() + 1);
        std::uniform_real_distribution<double> dist(0.0, 1.0);

#pragma omp for
        for (int64_t i = 0; i < n_boot; i++) {
            // refill the sample array with fresh draws
            for (int64_t j = 0; j < n; j++) {
                p_sample[j] = dist(rng);
            }
            // sort the allocated block for this bootstrap sample
            std::sort(p_sample, p_sample_end);
            dips[i] = dipstat <= diptst<false>(
                          p_sample,
                          n,
                          lo_hi.get(),
                          gcm.get(),
                          lcm.get(),
                          mn.get(),
                          mj.get(),
                          allow_zero,
                          debug);
        }
    }  // pragma parallel
    int64_t accu = 0;
    double p_val = static_cast<double>(
                       std::accumulate(dips.get(), dips.get() + n_boot, accu))
                   / n_boot;
    return p_val;
}  // diptest_pval_mt
#endif

namespace bindings {

void bind_diptest(py::module& m) {
    m.def(
        "diptest",
        &diptest::diptest,
        py::arg("x"),
        py::arg("allow_zero") = 1,
        py::arg("debug") = 0);
}

void bind_diptest_full(py::module& m) {
    m.def(
        "diptest_full",
        &diptest::diptest_full,
        py::arg("x"),
        py::arg("allow_zero") = 1,
        py::arg("debug") = 0);
}

void bind_diptest_pval(py::module& m) {
    m.def(
        "diptest_pval",
        &diptest::diptest_pval,
        py::arg("dipstat"),
        py::arg("n"),
        py::arg("n_boot") = 10000,
        py::arg("allow_zero") = 1,
        py::arg("debug") = 0,
        py::arg("seed") = 0,
        py::arg("stream") = 0);
}

#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
void bind_diptest_pval_mt(py::module& m) {
    m.def(
        "diptest_pval_mt",
        &diptest::diptest_pval_mt,
        py::arg("dipstat"),
        py::arg("n"),
        py::arg("n_boot") = 10000,
        py::arg("allow_zero") = 1,
        py::arg("debug") = 0,
        py::arg("seed") = 0,
        py::arg("n_threads") = 4);
}
#endif

}  // namespace bindings
}  // namespace diptest
