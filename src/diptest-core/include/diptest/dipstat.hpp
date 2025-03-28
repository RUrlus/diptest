/* dipstat.hpp -- header file for wrapper the diptest wrapper functions
 * Copyright 2022-2023 R. Urlus
 */
#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>  // unique_ptr

#include <diptest/common.hpp>
#include <diptest/diptest.hpp>

namespace py = pybind11;

namespace diptest {

namespace details {
double diptest(
    const double* x_ptr, int_vt N, int allow_zero = 1, int debug = 0
);
}  // namespace details

double diptest(const py::array_t<double>& x, int allow_zero, int debug);
py::dict diptest_full(const py::array_t<double>& x, int allow_zero, int debug);

namespace bindings {

void bind_diptest(py::module& m);
void bind_diptest_full(py::module& m);

}  // namespace bindings
}  // namespace diptest
