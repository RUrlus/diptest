/* wrapper.cpp -- implementation of wrapper around diptst from diptest.c
 * Copyright 2022 R. Urlus
 */
#include <cmath>   // log
#include <cstdint>
#include <memory>

#include <diptest/bootstrap.hpp>

namespace py = pybind11;

namespace diptest {

namespace details {

template <typename RNG>
inline double uniform_rvs(RNG& rng) {
    constexpr double inv_2p53 = 0x1.0p-53;  // 2^-53
    std::uint64_t r = rng();
    std::uint64_t k = r >> 11;
    return (static_cast<double>(k) + 1.0) * inv_2p53;
}

template <typename RNG>
inline double exp1(RNG& rng) {
    return -std::log(uniform_rvs(rng));
}

/*
 * Generate n sorted Uniform(0,1) samples using exponential spacings.
 *
 * This method exploits the fact that normalized partial sums of i.i.d.
 * Exponential(1) random variables have the same distribution as uniform
 * order statistics. This avoids O(n log n) sorting of uniform samples.
 *
 * Reference: Devroye, L. (1986). Non-Uniform Random Variate Generation.
 *            Springer-Verlag, Chapter 2, Section 2.2.
 */
template <typename RNG>
inline void sorted_uniform_rvs(double* out, int_vt n, RNG& rng) {
    double cumsum = 0.0;

    for (int_vt i = 0; i < n; ++i) {
        cumsum += exp1(rng);
        out[i] = cumsum;
    }
    cumsum += exp1(rng);

    const double inv_total = 1.0 / cumsum;
    for (int_vt i = 0; i < n; ++i) {
        out[i] *= inv_total;
    }
}

}  // namespace details

double diptest_pval(
    const double dipstat,
    const int64_t n,
    const int64_t n_boot,
    int allow_zero,
    int debug,
    uint64_t seed,
    uint64_t stream
) {
    details::pcg64_dxsm rng;
    if (seed == 0) {
        details::pcg_seed_seq seed_source;
        rng.seed(seed_source);
    } else if (stream != 0) {
        rng.seed(seed, stream);
    } else {
        rng.seed(seed);
    }

    std::array<int_vt, 5> lo_hi = {0, 0, 0, 0, 0};
    std::unique_ptr<int_vt[]> gcm(new int_vt[n]);
    std::unique_ptr<int_vt[]> lcm(new int_vt[n]);
    std::unique_ptr<int_vt[]> mn(new int_vt[n]);
    std::unique_ptr<int_vt[]> mj(new int_vt[n]);
    std::unique_ptr<double[]> sample(new double[n]);

    double* r_sample = sample.get();

    int64_t dip_cnt = 0;
    for (int64_t i = 0; i < n_boot; i++) {
        details::sorted_uniform_rvs(r_sample, n, rng);
        double dip = diptst<false>(
            r_sample,
            n,
            lo_hi.data(),
            gcm.get(),
            lcm.get(),
            mn.get(),
            mj.get(),
            allow_zero,
            debug
        );
        dip_cnt += dipstat <= dip;
    }
    return static_cast<double>(dip_cnt) / n_boot;
}  // diptest_pval

#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
double diptest_pval_mt(
    const double dipstat,
    const int64_t n,
    const int64_t n_boot,
    int allow_zero,
    int debug,
    uint64_t seed,
    size_t n_threads
) {
    details::pcg64_dxsm global_rng;
    if (seed == 0) {
        details::pcg_seed_seq seed_source;
        global_rng.seed(seed_source);
    } else {
        global_rng.seed(seed);
    }

    int64_t dip_cnt = 0;
#pragma omp parallel num_threads(n_threads) shared(global_rng)
    {
        std::array<int_vt, 5> lo_hi = {0, 0, 0, 0, 0};
        std::unique_ptr<int_vt[]> gcm(new int_vt[n]);
        std::unique_ptr<int_vt[]> lcm(new int_vt[n]);
        std::unique_ptr<int_vt[]> mn(new int_vt[n]);
        std::unique_ptr<int_vt[]> mj(new int_vt[n]);
        std::unique_ptr<double[]> sample(new double[n]);

        double* p_sample = sample.get();

        // PCG family has different streams which are, in theory, independent of
        // each other. Hence, we can use the same seed and a different stream to
        // draw independent samples from each thread without having to allocate
        // the whole block
        details::pcg64_dxsm rng = global_rng;
        rng.set_stream(omp_get_thread_num() + 1);

#pragma omp for reduction(+ : dip_cnt)
        for (int64_t i = 0; i < n_boot; i++) {
            details::sorted_uniform_rvs(p_sample, n, rng);
            dip_cnt += dipstat <= diptst<false>(
                           p_sample,
                           n,
                           lo_hi.data(),
                           gcm.get(),
                           lcm.get(),
                           mn.get(),
                           mj.get(),
                           allow_zero,
                           debug
                       );
        }
    }  // pragma parallel
    return static_cast<double>(dip_cnt) / n_boot;
}  // diptest_pval_mt
#endif

namespace bindings {

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
        py::arg("stream") = 0
    );
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
        py::arg("n_threads") = 4
    );
}
#endif

}  // namespace bindings
}  // namespace diptest
