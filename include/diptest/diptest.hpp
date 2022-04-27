/**
 * @file diptest.hpp
 * @author Prodromos Kolyvakis (prokolyvakis@gmail.com)
 * @brief header file including classes and methods helpful for dip calculation
 * @version 0.1
 * @date 2022-04-19
 *
 * @copyright Copyright (c) 2022 Prodromos Kolyvakis
 *
 */
#ifndef INCLUDE_DIPTEST_DIPTEST_HPP_
#define INCLUDE_DIPTEST_DIPTEST_HPP_
#define UNUSED(expr)  \
    do {              \
        (void)(expr); \
    } while (0)

#include <cassert>   // for assert
#include <cmath>     // for isgreaterequal
#include <iomanip>   // for setw
#include <iostream>  // for cout
#include <iterator>  // for iterators
#include <vector>    // for vectors

/**
 * @brief Enumerates the distinct types of functions that will be calculated
 * during the dip calculation
 *
 */
enum ConvexEnvelopeType { MAJORANT, MINORANT };

/**
 * @brief Storage of the dip value and the index at which the value was reported.
 *
 * The class, for the ease of speed, does not perform any consistency checks,
 * such as enforcing val and idx to be non-negative.
 *
 * @param val the dip value
 * @param idx the index at which the dip value is reported
 */
class Dip {
 public:
    double val;
    int idx;

    // Constructors:

    Dip(double val, int idx) : val(val), idx(idx) {}

    /**
     * @brief Construct a new Dip object
     *
     * @param min_is_0 if true the dip val will be set to 0. else to 1.
     */
    explicit Dip(bool min_is_0) : idx(-1) { val = (min_is_0) ? 0. : 1.; }

    // Methods:

    /**
     * @brief An update operation that makes the two instances to have the same
     * member variables if the parameter's dip value is greater the current one
     *
     * @param value the value that the stored value will be compared against
     * @param index the index that the dip value is reported
     */
    void maybe_update(double value, int index);

    /**
     * @brief An update operation that makes the two instances to have the same
     * member variables if the parameter's dip value is greater the current one
     *
     * @param other a Dip instance whose member parameters will be copied
     *
     * @overload
     */
    void maybe_update(const Dip& other);
};

inline void Dip::maybe_update(double value, int index) {
    if (val < value) {
        val = value;
        idx = index;
    }
}

inline void Dip::maybe_update(const Dip& other) {
    if (val < other.val) {
        val = other.val;
        idx = other.idx;
    }
}

/**
 * @brief A structure storing the needed parameters for the gcm (or lcm) fit
 *
 * @param arr the sorted array over which either the gcm or lcm will be computed
 * @param size the arr size
 * @param optimum an array storing either the gcm or lcm depending on the type
 * @param indices an array storing the indices needed for the optimum fit
 * @param rel_length the relevant length of either the gcm or lcm fit
 * @param x a counter for the convex majorant (or minorant)
 * @param y a counter for the convex majorant (or minorant)
 * @param type the type of the fit, i.e., either gcm or lcm fit
 */
class ConvexEnvelope {
 public:
    const double* arr;
    int *optimum, *indices;
    const int size;
    const ConvexEnvelopeType type;
    int rel_length = -1, x = -1, y = -1;

    // Constructors:

    ConvexEnvelope(const double* arr, int* optimum, int* indices, int size, ConvexEnvelopeType type)
        : arr(arr), optimum(optimum), indices(indices), size(size), type(type) {}

    // Methods:

    /**
     * @brief Establish the indices that are necessary for the convex minorant
     * (or majorant) fit
     *
     */
    void compute_indices();

    /**
     * @brief Computes the dip in the convex minorant (or majorant)
     *
     * @return the dip value and the index over which the dip was reported in
     * the array
     */
    Dip compute_dip();
};

/**
 * @brief Computes the greatest distance between the gcm and lcm
 *
 * @param gcm the current greatest convex minorant fit
 * @param lcm the current least convex majorant fit
 * @param debug the debug level. It is only relevant if it was compiled with
 * `DDIPTEST_ENABLE_DEBUG=ON`
 * @return the maximum distance
 */
double max_distance(ConvexEnvelope& gcm, ConvexEnvelope& lcm, int debug);

/**
 * @brief Calculates the dip for an ordered vector X using the greatest convex
 * minorant and the least concave majorant
 *
 * @param[in] x the array over which either the gcm or lcm will be computed
 * @param[in] n the size of the array
 * @param[out] lo_hi an array of size 4 that is used to return the lower and the
 * upper end of the model interval, and the relative lengths of gcm and lcm
 * @param[out] ifault an error integer. A value of 1 indicates that n is non-
 * positive. A value of 2 indicates that the array x was not sorted
 * @param gcm[out] the greatest convex minorant
 * @param lcm[out] the lowest convex majorant
 * @param mn[out] the greatest convex minorant's indices
 * @param mj[out] the greatest convex majorant's indices
 * @param min_is_0[in] a value indicating which is the dip's minimum value. If
 * set to 1, the minimum dip value can be 1., otherwise 0.
 * @param debug[in] the debug level. It is only relevant if it was compiled with
 * `DDIPTEST_ENABLE_DEBUG=ON`
 * @return the dip value
 */
double diptst(
    const double x[],
    const int n,
    int* lo_hi,
    int* ifault,
    int* gcm,
    int* lcm,
    int* mn,
    int* mj,
    const int min_is_0,
    const int debug);

#endif  // INCLUDE_DIPTEST_DIPTEST_HPP_
