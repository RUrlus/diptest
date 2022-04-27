/**
 * @file diptest.cpp
 * @author Martin Maechler (maechler@stat.math.ethz.ch)
 * @author Yong Lu (lyongu+@cs.cmu.edu)
 * @author Ralph Urlus (rurlus.dev@gmail.com)
 * @author Prodromos Kolyvakis (prokolyvakis@gmail.com)
 * @brief It performs the dip calculation for an ordered vector X using the
 * greatest convex minorant and the least concave majorant, skipping through
 * the data using the change points of these distributions.
 *
 * It returns the dip statistic 'DIP' and the modal interval (XL, XU).
 *                ===                        ======
 *
 * For further information regarding the Dip statistic please refer to the:
 *          ALGORITHM AS 217 APPL. STATIST. (1985) VOL.34, NO.3
 *
 * and, specifically, to the following articles:
 *
 * @article{HarP85,
 *   author = {P. M. Hartigan},
 *   title = {Computation of the Dip Statistic to Test for Unimodality},
 *   year = 1985,
 *   journal = {Applied Statistics},
 *   pages = {320--325},
 *    volume = 34 }
 *
 * @article{HarJH85,
 *   author = {J. A. Hartigan and P. M. Hartigan},
 *   title = {The Dip Test of Unimodality},
 *   year = 1985,
 *   journal = {Ann. of Statistics},
 *   pages = {70--84},
 *   volume = 13 }
 *
 *                              Historical Notes:
 *                              =================
 *
 * 1. The code was based on a Fortran version (dip.f) that was translated by
 *    f2c (version of 22 July 1992  22:54:52).
 * 2. The code was pretty-edited and extended (the debug argument) by:
 *    Martin Maechler <maechler@stat.math.ethz.ch>
 *     ETH Seminar fuer Statistik
 *     8092 Zurich    SWITZERLAND
 * 3. The overhead caused by debug branching statements was removed by placing
 *    them under a compile-time definition by:
 *    Ralph Urlus <rurlus.dev@gmail.com>
 * 4. A major code refactoring and transition from C to C++ was performed by:
 *    Prodromos Kolyvakis <prokolyvakis@gmail.com>
 *
 *                                Bug Fixes:
 *                                ==========
 *
 * 1. July 30 1994 : For unimodal data, gave "infinite loop"  (end of code)
 * 2. Oct  31 2003 : Yong Lu <lyongu+@cs.cmu.edu> : ")" typo in Fortran
 *                   gave wrong result (larger dip than possible) in some cases
 * 3. Apr  15 2022 : Prodromos Kolyvakis <prokolyvakis@gmail.com> : ">=" for
 *                   double comparison changed to the safer isgreaterequal
 *
 * @version 0.1
 * @date 2022-04-19
 *
 */
#include <diptest/diptest.hpp>

inline void ConvexEnvelope::compute_indices() {
    const int offset = (type == MINORANT) ? +1 : -1;
    const int start = (type == MINORANT) ? 1 : size;
    const int end = size + 1 - start;

    indices[start] = start;

    for (int i = start + offset; offset * (end - i) >= 0; i += offset) {
        indices[i] = i - offset;

        while (true) {
            int ind_at_i = indices[i];
            int ind_at_i_iter = indices[ind_at_i];

            /**
             * We compare the rate of change of arr, i.e., (arr[x]-arr[y])/(x-y),
             * at the indices:
             *      a. (i, ind_at_i)
             *      b. (ind_at_i, ind_at_i_iter)
             */
            bool rate_change_flag = (arr[i] - arr[ind_at_i]) * (ind_at_i - ind_at_i_iter)
                                    < (arr[ind_at_i] - arr[ind_at_i_iter]) * (i - ind_at_i);

            if (ind_at_i == start || rate_change_flag)
                break;
            indices[i] = ind_at_i_iter;
        }
    }
}

inline Dip ConvexEnvelope::compute_dip() {
    const int offset = (type == MINORANT) ? 0 : 1;
    const int sign = 1 + -2 * offset;

    Dip ret_dip(0., -1);
    Dip tmp_dip(1., -1);

    for (int j = x; j < rel_length; ++j) {
        int j_start = optimum[j + 1 - offset], j_end = optimum[j + offset];

        if (j_end - j_start > 1 && arr[j_end] != arr[j_start]) {
            double C = (j_end - j_start) / (arr[j_end] - arr[j_start]);

            for (int jj = j_start; jj <= j_end; ++jj) {
                double d = sign * ((jj - j_start + sign) - (arr[jj] - arr[j_start]) * C);

                tmp_dip.maybe_update(d, jj);
            }
        }

        ret_dip.maybe_update(tmp_dip);
        tmp_dip.val = 1.;
        tmp_dip.idx = -1;
    }

    return ret_dip;
}

double max_distance(ConvexEnvelope& gcm, ConvexEnvelope& lcm, int debug) {
#ifndef DDIPTEST_ENABLE_DEBUG
    UNUSED(debug);
#endif  // DDIPTEST_ENABLE_DEBUG

    assert(gcm.type != lcm.type && gcm.type == MINORANT);

    const double* arr = gcm.arr;
    long double ret_d = 0.;

    do {
        int gcm_y = gcm.optimum[gcm.y], lcm_y = lcm.optimum[lcm.y];
        int is_maj = gcm_y > lcm_y;
        int i = is_maj * gcm_y + (1-is_maj) * lcm_y;
        int j = is_maj * lcm_y + (1-is_maj) * gcm_y;
        int i1 = is_maj * (gcm.optimum[gcm.y + 1]) + (1-is_maj) * (lcm.optimum[lcm.y - 1]); 
        int sign = 2 * is_maj - 1;

        long double dx = sign * ((j - i1 + sign) - ((long double)arr[j] - arr[i1]) * (i - i1) / (arr[i] - arr[i1]));
        gcm.y -= (1 - is_maj);
        lcm.y += is_maj;

        if (std::isgreaterequal(dx, ret_d)) {
            ret_d = dx;
            gcm.x = gcm.y + 1;
            lcm.x = lcm.y - is_maj;
#if defined(DIPTEST_DEBUG)
            if (debug >= 2) {
                cout << ((is_maj) ? "G" : "L") << "(" << (gcm.x) << ", " << (lcm.x) << ")";
            }
#endif  // DIPTEST_DEBUG
        }
        if (gcm.y < 1)
            gcm.y = 1;
        if (lcm.y > lcm.rel_length)
            lcm.y = lcm.rel_length;

#if defined(DIPTEST_DEBUG)
        if (debug) {
            if (debug >= 2) {
                cout << " --> (gcm.y, lcm.y) = (" << (gcm.y) << ", " << (lcm.y) << ")" << endl;
            } else {
                cout << ".";
            }
        }
#endif  // DIPTEST_DEBUG
    } while (gcm.optimum[gcm.y] != lcm.optimum[lcm.y]);

    return ret_d;
}

// Subroutine
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
    const int debug) {
/**
 * `low` contains the index of the current estimate  of the lower end.
 * of the modal interval, `high` contains the index for the upper end.
 * It can be used as: double xl = x[low], xu = x[high];
 */
#define low lo_hi[0]
#define high lo_hi[1]

#ifndef DIPTEST_DEBUG
    UNUSED(debug);
#endif  // DIPTEST_DEBUG
    long double d = 0.;  // TODO: check if this makes 32/64-bit differences go
    double dip = (min_is_0) ? 0. : 1.;
    Dip dip_l(min_is_0), dip_u(min_is_0), tmp_dip(min_is_0);
    int i;
    bool flag;

    /**
     *  Parameter adjustments, so that array referencing starts at 1,
     *  i.e., x[1]..x[n]
     */
    --mj;
    --mn;
    --lcm;
    --gcm;
    --x;

    ConvexEnvelope gcm_obj(x, gcm, mn, n, MINORANT);
    ConvexEnvelope lcm_obj(x, lcm, mj, n, MAJORANT);

    // Perform the two consistency checks:

    // A. non-positive check:
    if (n <= 0) {
        *ifault = 1;
        return 0.0;
    }

    // B. non-sorted array check:
    for (i = 2; i <= n; ++i)
        if (x[i] < x[i - 1]) {
            *ifault = 2;
            return 0.0;
        }

    low = 1;
    high = n;

    if (n < 2 || x[n] == x[1])
        goto L_END;

#if defined(DIPTEST_DEBUG)
    if (debug)
        cout << "'dip': START: (N = " << n << ")"
             << " and 2N*dip = " << dip << "." << endl;
#endif  // DIPTEST_DEBUG

    /**
     * Establish the indices   mn[1..n]  over which combination is necessary
     * for the convex MINORANT (GCM) fit.
     */
    gcm_obj.compute_indices();

    /**
     * Establish the indices   mj[1..n]  over which combination is necessary
     * for the concave MAJORANT (LCM) fit.
     */
    lcm_obj.compute_indices();

    // ------------------------- Start the cycling. -------------------------
    do {
        // Collect the change points for the GCM from `high` to `low`.
        gcm[1] = high;
        for (i = 1; gcm[i] > low; i++)
            gcm[i + 1] = mn[gcm[i]];
        gcm_obj.x = gcm_obj.rel_length = i;  //< relevant_length(GCM)
        gcm_obj.y = gcm_obj.x - 1;

        // Collect the change points for the LCM from `high` to `low`.
        lcm[1] = low;
        for (i = 1; lcm[i] < high; i++)
            lcm[i + 1] = mj[lcm[i]];
        lcm_obj.x = lcm_obj.rel_length = i;  //< relevant_length(GCM)
        lcm_obj.y = 2;

/**
 * - `l_gcm` is defined as: relevant_length(GCM).
 * - `l_lcm` is defined as: relevant_length(LCM)
 */
#define l_gcm (gcm_obj.rel_length)
#define l_lcm (lcm_obj.rel_length)

#if defined(DIPTEST_DEBUG)
        if (debug) {
            cout << "'dip': LOOP-BEGIN: 2n*D = " << dip << " and [low, high] = [" << setw(3) << low << ", " << setw(3)
                 << high << "]";
            if (debug >= 3) {
                // Print the GCM:
                cout << " :" << endl << " gcm[1:" << l_gcm << "] = ";
                for (i = 1; i < l_gcm; i++)
                    cout << gcm[i] << ", ";
                cout << gcm[l_gcm] << endl;
                // Print the LCM:
                cout << " lcm[1:" << l_lcm << "] = ";
                for (i = 1; i < l_lcm; i++)
                    cout << lcm[i] << ", ";
                cout << lcm[l_lcm] << endl;
            } else {  // debug <= 2
                cout << "; (l_lcm, l_gcm) = (" << setw(2) << l_lcm << ", " << setw(2) << l_gcm << ")" << endl;
            }
        }
#endif  // DIPTEST_DEBUG
        /**
         * Find the largest distance greater than 'DIP' between
         * the GCM and the LCM from LOW to HIGH.
         */
        if (l_gcm != 2 || l_lcm != 2) {
#if defined(DIPTEST_DEBUG)
            if (debug) {
                cout << "'dip': CYCLE-BEGIN: while(gcm[gcm_obj.y] "
                     << "!= lcm[lcm_obj.y])";
                if (debug >= 2) {
                    cout << endl;
                } else {
                    cout << " ";
                }
            }
#endif  // DIPTEST_DEBUG

            d = max_distance(gcm_obj, lcm_obj, debug);

#if defined(DIPTEST_DEBUG)
            if (debug && debug < 2)
                cout << endl;
#endif  // DIPTEST_DEBUG
        } else {
            d = (min_is_0) ? 0. : 1.;

#if defined(DIPTEST_DEBUG)
            if (debug)
                cout << "'dip': NO-CYCLE: (l_lcm, l_gcm) = (" << setw(2) << l_lcm << ", " << setw(2) << l_gcm
                     << ") ==> d := " << d << endl;
#endif  // DIPTEST_DEBUG
        }

        if (d < dip)
            goto L_END;

        // Calculate the DIPs for the current LOW and HIGH
#if defined(DIPTEST_DEBUG)
        if (debug)
            cout << "'dip': MAIN-CALCULATION" << endl;
#endif  // DIPTEST_DEBUG

        // The DIP for the convex minorant.
        dip_l = gcm_obj.compute_dip();

        // The DIP for the concave majorant.
        dip_u = lcm_obj.compute_dip();

#if defined(DIPTEST_DEBUG)
        if (debug)
            cout << " (dip_l, dip_u) = (" << dip_l.val << ", " << dip_u.val << ")";
#endif  // DIPTEST_DEBUG

        // Determine the current maximum.
        if (dip_l.val < dip_u.val) {
            tmp_dip = dip_u;
        } else {
            tmp_dip = dip_l;
        }
        if (dip < tmp_dip.val) {
            dip = tmp_dip.val;
#if defined(DIPTEST_DEBUG)
            if (debug)
                cout << " --> new larger dip " << (tmp_dip.val) << " ( at index := " << (tmp_dip.idx) << " )" << endl;
#endif  // DIPTEST_DEBUG
        }

        flag = (low == gcm[gcm_obj.x] && high == lcm[lcm_obj.x]);
        low = gcm[gcm_obj.x];
        high = lcm[lcm_obj.x];
    } while (!flag);
    // ------------------------- end of the cycling. -------------------------

#if defined(DIPTEST_DEBUG)
    if (debug)
        cout << "'dip': LOOP-END: No improvement found neither in low := " << low << " nor in high := " << high << endl;
#endif  // DIPTEST_DEBUG

L_END:
    /**
     * M. Maechler -- speedup: Work with (2n * dip) everywhere but the very end!
     * It saves many divisions by n!
     */
    dip /= (2 * n);
    lo_hi[2] = l_gcm;
    lo_hi[3] = l_lcm;
    return dip;
}  // diptst
#undef low
#undef high
#undef l_gcm
#undef l_lcm
