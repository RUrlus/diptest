/* pcg64.hpp -- Typedefs of DXSM variant of PCG64
 * Copyright 2022 R. Urlus */
#ifndef INCLUDE_DIPTEST_PCG64_HPP_
#define INCLUDE_DIPTEST_PCG64_HPP_

#include <random>  // random_device

#include <pcg_extras.hpp>
#include <pcg_random.hpp>

namespace diptest {

typedef pcg_engines::setseq_dxsm_128_64 pcg64_dxsm;
typedef pcg_extras::seed_seq_from<std::random_device> pcg_seed_seq;

}  // namespace diptest

#endif  // INCLUDE_DIPTEST_PCG64_HPP_
