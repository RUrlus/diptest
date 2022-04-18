#ifndef INCLUDE_DIPTEST_DIP_HPP_
#define INCLUDE_DIPTEST_DIP_HPP_
#define UNUSED(expr)  \
    do                \
    {                 \
        (void)(expr); \
    } while (0)

#include <iterator> // for iterators
#include <vector>   // for vectors
#include <assert.h> // for assert
#include <iostream> // for cout
#include <iomanip>  // for setw
#include <math.h>   // isgreaterequal

using namespace std;

enum ConvexEnvelopeType
{
    MAJORANT,
    MINORANT
};

class DipValue
{
public:
    double val;
    int idx;

    // Constructor:
    DipValue(double val, int idx) : val(val), idx(idx){};
    DipValue(bool min_is_0) : idx(-1) { val = (min_is_0) ? 0. : 1.; };

    // Methods:
    void update(double value, int index)
    {
        val = value;
        idx = index;
    }

    void update(DipValue &other)
    {
        val = other.val;
        idx = other.idx;
    }

    void maybe_update(double value, int index)
    {
        if (val < value)
        {
            val = value;
            idx = index;
        }
    }

    void maybe_update(DipValue &other)
    {
        maybe_update(other.val, other.idx); 
    }
};

class ConvexEnvelope
{

public:
    const double *arr;
    const int size;
    int *optimum, *indices;
    int rel_length = -1, x = -1, y = -1;
    const ConvexEnvelopeType type;

    ConvexEnvelope(
        const double *arr,
        int *optimum,
        int *indices,
        int size,
        ConvexEnvelopeType type) : arr(arr), optimum(optimum), indices(indices), size(size), type(type){};

    /* Establish the indices   mn[1..n]  over which combination is necessary
       for the convex MINORANT (GCM) fit.
     */
    void compute_indices();

    DipValue compute_dip();
};

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
    const int debug);

#endif // INCLUDE_DIPTEST_DIP_HPP_