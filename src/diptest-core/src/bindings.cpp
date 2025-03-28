/* bindings.cpp -- Python bindings for diptest
 * Copyright 2022 Ralph Urlus
 */
#include <pybind11/pybind11.h>
#include <diptest/bootstrap.hpp>
#include <diptest/dipstat.hpp>
#if defined(DIPTEST_BUILD_CPP_TESTS)
#include <diptest/test_pcg.hpp>
#endif  // DIPTEST_BUILD_CPP_TESTS

namespace diptest {
namespace bindings {

PYBIND11_MODULE(_diptest_core, m) {
    bind_diptest(m);
    bind_diptest_full(m);
    bind_diptest_pval(m);
#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
    bind_diptest_pval_mt(m);
#endif

#if defined(DIPTEST_BUILD_CPP_TESTS)
    bind_pcg_seed_test(m);
    bind_pcg_set_stream_test(m);
#if defined(DIPTEST_HAS_OPENMP_SUPPORT)
    bind_pcg_mt_stream_test(m);
#endif  // DIPTEST_HAS_OPENMP_SUPPORT
#endif  // DIPTEST_BUILD_CPP_TESTS

#ifndef OS_WIN
#ifdef DIPTEST_VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(DIPTEST_VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
#endif

#ifdef DIPTEST_HAS_OPENMP_SUPPORT
    m.attr("_has_openmp_support") = true;
#else
    m.attr("_has_openmp_support") = false;
#endif  // MMU_HAS_OPENMP_SUPPORT
}

}  // namespace bindings
}  // namespace diptest
