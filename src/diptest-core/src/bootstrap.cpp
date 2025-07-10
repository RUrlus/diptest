/* wrapper.cpp -- implementation of wrapper around diptst from diptest.c
 * Copyright 2022 R. Urlus
 */
#include <algorithm>  // sort
#include <memory>
#include <numeric>  // accumulate
#include <random>   // uniform_real_distribution

#include <diptest/bootstrap.hpp>

namespace py = pybind11;

namespace diptest {

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
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    std::array<int_vt, 5> lo_hi = {0, 0, 0, 0, 0};
    std::unique_ptr<int_vt[]> gcm(new int_vt[n]);
    std::unique_ptr<int_vt[]> lcm(new int_vt[n]);
    std::unique_ptr<int_vt[]> mn(new int_vt[n]);
    std::unique_ptr<int_vt[]> mj(new int_vt[n]);
    std::unique_ptr<int_vt[]> dips(new int_vt[n_boot]);
    std::unique_ptr<double[]> sample(new double[n]);

    double* r_sample = sample.get();
    double* sample_end = r_sample + n;

    for (int64_t i = 0; i < n_boot; i++) {
        for (int64_t j = 0; j < n; j++) {
            r_sample[j] = dist(rng);
        }
        std::sort(r_sample, sample_end);
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
        dips[i] = dipstat <= dip;
    }
    int64_t accu = 0;
    double p_val = static_cast<double>(
                       std::accumulate(dips.get(), dips.get() + n_boot, accu)
                   )
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
    size_t n_threads
) {
    std::unique_ptr<bool[]> dips(new bool[n_boot]);
    details::pcg64_dxsm global_rng;
    if (seed == 0) {
        details::pcg_seed_seq seed_source;
        global_rng.seed(seed_source);
    } else {
        global_rng.seed(seed);
    }

#pragma omp parallel num_threads(n_threads) shared(dips, global_rng)
    {
        std::array<int_vt, 5> lo_hi = {0, 0, 0, 0, 0};
        std::unique_ptr<int_vt[]> gcm(new int_vt[n]);
        std::unique_ptr<int_vt[]> lcm(new int_vt[n]);
        std::unique_ptr<int_vt[]> mn(new int_vt[n]);
        std::unique_ptr<int_vt[]> mj(new int_vt[n]);
        std::unique_ptr<double[]> sample(new double[n]);

        double* p_sample = sample.get();
        double* p_sample_end = p_sample + n;

        // PCG family has different streams which are, in theory, independent of
        // each other. Hence, we can use the same seed and a different stream to
        // draw independent samples from each thread without having to allocate
        // the whole block
        details::pcg64_dxsm rng = global_rng;
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
    int64_t accu = 0;
    double p_val = static_cast<double>(
                       std::accumulate(dips.get(), dips.get() + n_boot, accu)
                   )
                   / n_boot;
    return p_val;
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
