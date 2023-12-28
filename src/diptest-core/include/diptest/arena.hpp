#pragma once

#include <array>
#include <diptest/common.hpp>
#include <memory>

namespace diptest {
namespace details {

struct DiptestArena {
    pcg64_dxsm rng;
    std::uniform_real_distribution<double> dist;
    std::array<int, 4> lo_hi = {0, 0, 0, 0};
    std::unique_ptr<int[]> gcm;
    std::unique_ptr<int[]> lcm;
    std::unique_ptr<int[]> mn;
    std::unique_ptr<int[]> mj;
    std::unique_ptr<int[]> dips;
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
        uint64_t stream)
        : dist{std::uniform_real_distribution<double>(0.0, 1.0)},
          gcm{std::unique_ptr<int[]>(new int[n])},
          lcm{std::unique_ptr<int[]>(new int[n])},
          mn{std::unique_ptr<int[]>(new int[n])},
          mj{std::unique_ptr<int[]>(new int[n])},
          dips{std::unique_ptr<int[]>(new int[n_boot])},
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
