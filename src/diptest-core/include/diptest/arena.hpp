#pragma once

#include <array>
#include <memory>

#include <diptest/common.hpp>

namespace diptest {
namespace details {

struct DiptestArena {
    pcg64_dxsm rng;
    std::uniform_real_distribution<double> dist;
    std::array<int_vt, 4> lo_hi = {0, 0, 0, 0};
    std::unique_ptr<int_vt[]> gcm;
    std::unique_ptr<int_vt[]> lcm;
    std::unique_ptr<int_vt[]> mn;
    std::unique_ptr<int_vt[]> mj;
    std::unique_ptr<int_vt[]> dips;
    std::unique_ptr<double[]> sample;
    double* sample_begin;
    const int64_t n_boot;
    int allow_zero;
    int debug;

    DiptestArena(
        const int64_t n,
        const int64_t n_boot,
        int allow_zero,
        int debug,
        uint64_t seed,
        uint64_t stream
    )
        : dist{std::uniform_real_distribution<double>(0.0, 1.0)},
          gcm{std::unique_ptr<int_vt[]>(new int_vt[n])},
          lcm{std::unique_ptr<int_vt[]>(new int_vt[n])},
          mn{std::unique_ptr<int_vt[]>(new int_vt[n])},
          mj{std::unique_ptr<int_vt[]>(new int_vt[n])},
          dips{std::unique_ptr<int_vt[]>(new int_vt[n_boot])},
          sample{std::unique_ptr<double[]>(new double[n])},
          sample_begin{sample.get()},
          n_boot{n_boot},
          allow_zero{allow_zero},
          debug{debug} {
        if (seed == 0) {
            pcg_seed_seq seed_source;
            rng.seed(seed_source);
        } else if (stream != 0) {
            rng.seed(seed, stream);
        } else {
            rng.seed(seed);
        }
    }

};  // DiptestArena
//
}  // namespace details
}  // namespace diptest
