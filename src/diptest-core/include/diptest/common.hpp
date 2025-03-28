/* wrapper.hpp -- header file for wrapper around diptest implementation
 * Copyright 2022 R. Urlus
 */
#pragma once

#if defined(_WIN32) || defined(_WIN64) || defined(WIN32) \
    || defined(__CYGWIN__) || defined(__MINGW32__) || defined(__BORLANDC__)
#define OS_WIN
#endif

// handle error C2059: syntax error: ';'  on windows for this Macro
#ifndef OS_WIN
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
#endif

// Fix for lack of ssize_t on Windows for CPython3.10
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4127 \
)  // warning C4127: Conditional expression is constant
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

#include <cstdint>
#include <random>  // random_device

#include <pcg_extras.hpp>
#include <pcg_random.hpp>

#if defined(DIPTEST_64BIT_INDEX)
using int_vt = int64_t;
#else
using int_vt = int32_t;
#endif

namespace diptest {
namespace details {

typedef pcg_engines::cm_setseq_dxsm_128_64 pcg64_dxsm;
typedef pcg_extras::seed_seq_from<std::random_device> pcg_seed_seq;

}  // namespace details
}  // namespace diptest
