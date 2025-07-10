/* wrapper.cpp -- implementation of wrapper around diptst from diptest.c
 * Copyright 2022 R. Urlus
 */
#include <diptest/dipstat.hpp>

namespace py = pybind11;

namespace diptest {
namespace details {

inline double diptest(
    const double* x_ptr, int_vt N, int allow_zero, int debug
) {
    std::array<int_vt, 5> lo_hi = {0, 0, 0, 0, 0};
    std::unique_ptr<int_vt[]> gcm(new int_vt[N]);
    std::unique_ptr<int_vt[]> lcm(new int_vt[N]);
    std::unique_ptr<int_vt[]> mn(new int_vt[N]);
    std::unique_ptr<int_vt[]> mj(new int_vt[N]);

    double dip = diptst<true>(
        x_ptr,
        N,
        lo_hi.data(),
        gcm.get(),
        lcm.get(),
        mn.get(),
        mj.get(),
        allow_zero,
        debug
    );
    return dip;
}  // diptest

}  // namespace details

double diptest(const py::array_t<double>& x, int allow_zero, int debug) {
    return details::diptest(x.data(), x.size(), allow_zero, debug);
}  // diptest

py::dict diptest_full(const py::array_t<double>& x, int allow_zero, int debug) {
    const double* x_ptr = x.data();
    int_vt N = x.size();
    std::array<int_vt, 5> lo_hi = {0, 0, 0, 0, 0};

    auto gcm = py::array_t<int_vt>(N);
    auto lcm = py::array_t<int_vt>(N);
    int_vt* gcm_ptr = gcm.mutable_data();
    int_vt* lcm_ptr = lcm.mutable_data();

    std::unique_ptr<int_vt[]> mn(new int_vt[N]);
    std::unique_ptr<int_vt[]> mj(new int_vt[N]);
    int_vt* mn_ptr = mn.get();
    int_vt* mj_ptr = mj.get();

    double dip = diptst<true>(
        x_ptr,
        N,
        lo_hi.data(),
        gcm_ptr,
        lcm_ptr,
        mn_ptr,
        mj_ptr,
        allow_zero,
        debug
    );

    using namespace pybind11::literals;  // to bring in the `_a` literal NOLINT
    // NOTE diptst uses indexing starting from 1, so all the indexes returned
    // need to be corrected
    int_vt lo = lo_hi[0] - 1;
    int_vt hi = lo_hi[1] - 1;
    return py::dict(
        "dip"_a = dip,
        "lo"_a = lo,
        "hi"_a = hi,
        "dipidx"_a = lo_hi[4] - 1,
        "xl"_a = x.at(lo),
        "xu"_a = x.at(hi),
        "_gcm"_a = gcm,
        "_lcm"_a = lcm,
        "_lh_2"_a = lo_hi[2] - 1,
        "_lh_3"_a = lo_hi[3] - 1
    );
}  // diptest_full

namespace bindings {

void bind_diptest(py::module& m) {
    m.def(
        "diptest",
        &diptest::diptest,
        py::arg("x"),
        py::arg("allow_zero") = 1,
        py::arg("debug") = 0
    );
}

void bind_diptest_full(py::module& m) {
    m.def(
        "diptest_full",
        &diptest::diptest_full,
        py::arg("x"),
        py::arg("allow_zero") = 1,
        py::arg("debug") = 0
    );
}

}  // namespace bindings
}  // namespace diptest
