#ifndef INCLUDE_DIPTEST_DIP_HPP_
#define INCLUDE_DIPTEST_DIP_HPP_
#define UNUSED(expr) do { (void)(expr); } while (0)

#include <iterator> // for iterators
#include <vector>   // for vectors
#include <assert.h> // for assert
#include <iostream> // for cout
#include <iomanip>  // for setw

using namespace std;

void compute_convex_m_indices(const double *arr, int *ret_idx, const std::vector<int>& range);

void compute_dip(const double *arr, const int *convex_m, int rel_length, int idx, int offset, double *ret_dip, int *ret_dip_idx);

long double compute_largest_distance_greater_than_dip(const double *arr, const int *convex_maj, const int *convex_min, int *ig, int *ix, int *ih, int *iv, int l_lcm, int debug);

/* Subroutine */
double diptst(
    const double x[],
    const int n,
    int *lo_hi,
    int *ifault,
    int *gcm,
    int *lcm,
    int *mn,
    int *mj,
    const int min_is_0,
    const int debug
);

#endif  // INCLUDE_DIPTEST_DIP_HPP_