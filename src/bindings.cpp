/* bindings.cpp -- Python bindings for diptest
 * Copyright 2022 Ralph Urlus
 */
#include <pybind11/pybind11.h>
#include <diptest/wrapper.hpp>
#if defined(DIPTEST_BUILD_CPP_TESTS)
#include <diptest/test_pcg.hpp>
#endif // DIPTEST_BULD_CPP_TESTS

namespace diptest {
namespace bindings {

PYBIND11_MODULE(EXTENSION_MODULE_NAME, m) {
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
#endif // DIPTEST_HAS_OPENMP_SUPPORT
#endif // DIPTEST_BULD_CPP_TESTS

#ifndef OS_WIN
#ifdef DIPTEST_VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(DIPTEST_VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
#endif
}

}  // namespace bindings
}  // namespace diptest
