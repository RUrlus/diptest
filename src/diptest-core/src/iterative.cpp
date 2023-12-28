/* wrapper.cpp -- implementation of wrapper around diptst from diptest.c
 * Copyright 2022 R. Urlus
 */
#include <diptest/iterative.hpp>
#include <iostream>

namespace py = pybind11;

namespace diptest {
namespace details {

using ivT = std::pair<int, int>;

struct DiptestArena {
    pcg64_dxsm rng;
    std::uniform_real_distribution<double> dist;
    std::array<int, 4> lo_hi_arr = {0, 0, 0, 0};
    std::unique_ptr<int[]> arena;
    std::unique_ptr<double[]> sample;
    double* sample_begin;
    const int64_t n_boot;
    double dipstat;
    int* gcm;
    int* lcm;
    int* mn;
    int* mj;
    int* dips;
    int* lo_hi = lo_hi_arr.data();
    int allow_zero;
    int debug;

    DiptestArena(
        const int n,
        const int n_boot,
        int allow_zero,
        int debug,
        uint64_t seed,
        uint64_t stream)
        : dist{std::uniform_real_distribution<double>(0.0, 1.0)},
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
        auto arena = std::unique_ptr<int[]>(new int[4 * n + n_boot]);
        int* data = arena.get();
        gcm = data;
        data += n;
        lcm = data;
        data += n;
        mn = data;
        data += n;
        mj = data;
        data += n;
        dips = data;
        sample = std::unique_ptr<double[]>(new double[n]);
        sample_begin = sample.get();
    }

};  // DiptestArena

double bootstrap(const int n, DiptestArena& arena) {
    double* sample_end = arena.sample_begin + n;
    for (int i = 0; i < arena.n_boot; i++) {
        for (int j = 0; j < n; j++) {
            arena.sample_begin[j] = arena.dist(arena.rng);
        }
        std::sort(arena.sample_begin, sample_end);
        double dip = diptst<false>(
            arena.sample_begin,
            n,
            arena.lo_hi,
            arena.gcm,
            arena.lcm,
            arena.mn,
            arena.mj,
            arena.allow_zero,
            arena.debug);
        arena.dips[i] = arena.dipstat <= dip;
    }
    int64_t accu = 0;
    double p_val = static_cast<double>(std::accumulate(
                       arena.dips, arena.dips + arena.n_boot, accu))
                   / arena.n_boot;
    return p_val;
}

void diptest(const double* x_ptr, const int n, DiptestArena& arena) {
    arena.dipstat = diptst<false>(
        x_ptr,
        n,
        arena.lo_hi,
        arena.gcm,
        arena.lcm,
        arena.mn,
        arena.mj,
        arena.allow_zero,
        arena.debug);
}

struct IV {
    // the lower bound of the interval
    int lo;
    // the upper bound of the interval
    int hi;
    // the size of the interval
    int size;
    // the minimum value the lower bound can take
    const int min;
    // the maximum value the upper bound can take
    const int max;
    // the offset from the global zero idx
    const size_t offset;
    // the start of the global memory that contains at least the interval itself
    const double* const data;
    // begin points to the start of the memory of the interval, i.e. data +
    // offset + lo => begin
    const double* begin;

    IV(int lo, int hi, const IV& other)
        : lo{lo},
          hi{hi},
          size{1 + hi - lo},
          min{other.min},
          max{other.max},
          offset{other.offset + (other.lo - lo)},
          data{other.data},
          begin{data + offset + lo} {}
};

struct Interval {
    int lo;
    int hi;
    int iv_size;
    int lb;
    int ub;
    int b_size;
    // the start of the global memory that contains at least the interval itself
    const double* const data;
    // begin points to the start of the memory of the interval
    const double* begin;

    Interval(int lo, int hi, const double* data)
        : lo{lo},
          hi{hi},
          iv_size{1 + hi - lo},
          lb{lo},
          ub{hi},
          b_size{1 + ub - lb},
          data{data},
          begin{data + lo} {}

    Interval(int lo, int hi, const double* begin, const Interval& other)
        : lo{lo},
          hi{hi},
          iv_size{1 + hi - lo},
          lb{other.lb},
          ub{other.ub},
          b_size{other.b_size},
          data{other.data},
          begin{begin} {}

    Interval(int lo, int hi, const Interval& other)
        : lo{lo},
          hi{hi},
          iv_size{1 + hi - lo},
          lb{other.lb},
          ub{other.ub},
          b_size{other.b_size},
          data{other.data},
          begin{data + lo} {}
};

template <bool left>
Interval local_to_global(Interval iv) {
    if (left) {
        if (iv.begin[iv.hi] < 0.0) {
            return iv;
        }
    } else {
        if (iv.begin[iv.lo] < 1e-15) {
            return iv;
        }
    }
    const int shift = iv.b_size - 1;
    iv.lo -= shift;
    iv.hi -= shift;
    iv.begin = iv.data + iv.lo;
    return iv;
}

Interval obtain_full_modal_interval(
    Interval& interval,
    DiptestArena& arena,
    double* __restrict mirror_dest) {
    int64_t mirror_size = 2 * interval.iv_size - 1;
    const double* rhs_ptr = mirror_dest + mirror_size;
    mirror_array<double>(
        interval.b_size, mirror_size, interval.begin, mirror_dest);
    // diptest on left mirror
    double ldt = diptst<false>(
        mirror_dest,
        mirror_size,
        arena.lo_hi,
        arena.gcm,
        arena.lcm,
        arena.mn,
        arena.mj,
        arena.allow_zero,
        arena.debug);
    auto liv = Interval(arena.lo_hi[0], arena.lo_hi[1], mirror_dest, interval);
    double rdt = diptst<false>(
        rhs_ptr,
        mirror_size,
        arena.lo_hi,
        arena.gcm,
        arena.lcm,
        arena.mn,
        arena.mj,
        arena.allow_zero,
        arena.debug);
    if (ldt > rdt) {
        auto tmp = local_to_global<true>(liv);
        return tmp;
    }
    return local_to_global<false>(
        Interval(arena.lo_hi[0], arena.lo_hi[1], rhs_ptr, interval));
}

void create_intervals(std::queue<Interval>& queue, const Interval& interval) {
    if (interval.lo - interval.lb > 10) {
        queue.push(Interval(interval.lb, interval.lo, interval.data));
    }
    if (interval.ub - interval.hi > 10) {
        queue.push(Interval(interval.hi, interval.ub, interval.data));
    }
}

}  // namespace details
//
//
//
void print_queue(const std::queue<details::Interval>& src) {
    py::print("queue contains: \n");
    std::queue<details::Interval> dest = src;

    while (!dest.empty()) {
        auto iv = dest.front();
        py::print("[", iv.lo, ", ", iv.hi, "] (", iv.lb, ", ", iv.ub, ")\n");
        dest.pop();
    }
}

/* idip -- iterative diptest
 *
 * Find modal intervals by performing the diptest iteratively.
 * The function performs breath-first iteration over the range using a
 * std::queue
 *
 *
 */
py::array_t<double> iterdip(
    const py::array_t<double>& x,
    const double crit_val,
    const int n_boot,
    const int max_iter,
    const int allow_zero,
    const int debug,
    uint64_t seed,
    uint64_t stream) {
    const int mirror_size = (x.size() * 2) - 1;

    // allocate the memory used to store the mirrors of (subsets) x
    auto mirror_arena = std::unique_ptr<double[]>(new double[2 * mirror_size]);
    // allocate the memory needed for diptest making
    // sure allocate enough space for mirrors of the data
    double* mirror_dest = mirror_arena.get();
    auto arena = details::DiptestArena(
        mirror_size, n_boot, allow_zero, debug, seed, stream);
    const double* x_ptr = x.data();

    auto queue = std::queue<details::Interval>();
    auto intervals = std::vector<details::Interval>();
    intervals.reserve(30);

    // always start with the whole interval
    queue.push(details::Interval(0, x.size() - 1, x_ptr));
    print_queue(queue);

    int iter = 0;
    while ((queue.size() > 0) && (iter <= max_iter)) {
        details::Interval modal_iv = queue.front();
        details::diptest(x_ptr + modal_iv.lo, modal_iv.iv_size, arena);
        double pval = details::bootstrap(modal_iv.iv_size, arena);
        py::print("pval: ", pval);
        if (pval < crit_val) {
            // we reject the null hypothesis of unimodality
            // we constrict search space to the modal interval of the previous
            // iteration
            queue.push(
                details::Interval(arena.lo_hi[0], arena.lo_hi[1], modal_iv));

            py::print("pval <, added modal interval to queue");
        } else {
            // we extend the modal interval as it can be too small
            //
            details::Interval extended_miv
                = details::obtain_full_modal_interval(
                    modal_iv, arena, mirror_dest);

            py::print(
                "extended_miv: ",
                extended_miv.lo,
                " - ",
                extended_miv.hi,
                " bounds ",
                extended_miv.lb,
                " - ",
                extended_miv.ub,
                "\n");
            // add to the vector of final intervals
            intervals.push_back(extended_miv);
            if (iter > 0) {
                // create new intervals around the extended modal interval
                details::create_intervals(queue, extended_miv);
            }
            py::print("pval >");
        }
        iter++;
        queue.pop();
        py::print("post iter");
        print_queue(queue);
    }

    const int rN = intervals.size();
    auto result = py::array_t<double>(
        {py::ssize_t_cast(rN), static_cast<py::ssize_t>(2)});
    auto res_ptr = result.mutable_unchecked<2>();

    for (int i = 0; i < rN; ++i) {
        auto iv = intervals.at(i);
        res_ptr(i, 0) = x_ptr[iv.lo];
        // intervals are half open [iv.lo, iv.hi)
        res_ptr(i, 1) = x_ptr[iv.hi];
    }
    return result;
}

namespace bindings {

void bind_mirror_array(py::module& m) {
    m.def(
        "mirror_array",
        [](const py::array_t<int>& arr) {
            return diptest::mirror_array<int>(arr);
        },
        py::arg("arr").noconvert(),
        py::return_value_policy::take_ownership);
    m.def(
        "mirror_array",
        [](const py::array_t<int64_t>& arr) {
            return diptest::mirror_array<int64_t>(arr);
        },
        py::arg("arr").noconvert(),
        py::return_value_policy::take_ownership);
    m.def(
        "mirror_array",
        [](const py::array_t<float>& arr) {
            return diptest::mirror_array<float>(arr);
        },
        py::arg("arr").noconvert(),
        py::return_value_policy::take_ownership);
    m.def(
        "mirror_array",
        [](const py::array_t<double>& arr) {
            return diptest::mirror_array<double>(arr);
        },
        py::arg("arr"),
        py::return_value_policy::take_ownership);
}

void bind_iterdip(py::module& m) {
    m.def(
        "iterdip",
        &diptest::iterdip,
        py::arg("x"),
        py::arg("crit_val") = 0.1,
        py::arg("n_boot") = 10000,
        py::arg("max_iter") = 100,
        py::arg("allow_zero") = 1,
        py::arg("debug") = 0,
        py::arg("seed") = 0,
        py::arg("stream") = 0);
}

}  // namespace bindings
}  // namespace diptest
