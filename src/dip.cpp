/*   ALGORITHM AS 217 APPL. STATIST. (1985) VOL.34, NO.3

  @article{HarP85,
     author = {P. M. Hartigan},
     title = {Computation of the Dip Statistic to Test for Unimodality},
     year = 1985,
     journal = {Applied Statistics},
     pages = {320--325},
     volume = 34 }
  @article{HarJH85,
     author = {J. A. Hartigan and P. M. Hartigan},
     title = {The Dip Test of Unimodality},
     year = 1985,
     journal = {Ann. of Statistics},
     pages = {70--84},
     volume = 13 }

  Does the dip calculation for an ordered vector X using the
  greatest convex minorant and the least concave majorant, skipping
  through the data using the change points of these distributions.

  It returns the dip statistic 'DIP' and the modal interval (XL, XU).
                 ===                        ======

   dip.f -- translated by f2c (version of 22 July 1992  22:54:52).

   Pretty-Edited and extended (debug argument)
   by Martin Maechler <maechler@stat.math.ethz.ch>
      ETH Seminar fuer Statistik
      8092 Zurich    SWITZERLAND

---------------

   Two Bug Fixes:
       =========

   1)   July 30 1994 : For unimodal data, gave "infinite loop"  (end of code)
   2)   Oct  31 2003 : Yong Lu <lyongu+@cs.cmu.edu> : ")" typo in Fortran
                       gave wrong result (larger dip than possible) in some cases
   $Id: dip.c,v 1.26 2012/08/13 16:44:11 maechler Exp $
*/
#include <diptest/dip.hpp>

void compute_convex_m_indices(const double *arr, int *ret_idx, const std::vector<int>& range) {
    vector<int>::const_iterator curr_idx = range.begin();
    vector<int>::const_iterator next_idx = next(curr_idx, 1);
    
    int begin_v = *curr_idx;
    ret_idx[begin_v] = begin_v;

    for (; next_idx < range.end(); curr_idx++, next_idx++)  {
        ret_idx[*next_idx] = *curr_idx;
        
        while(1) {
            int next_v = ret_idx[*next_idx];
            int next_v_iter = ret_idx[next_v];

            /*
             * We compare the rate of change function of arr at the indices:
             *      a. (*next_idx, next_v)
             *      b. (next_v, next_v_iter)
             */
            bool rate_change_flag = 
                (arr[*next_idx]  - arr[next_v])  * (next_v - next_v_iter) <
                (arr[next_v] - arr[next_v_iter]) * (*next_idx - next_v);
            
            if (next_v_iter == begin_v || rate_change_flag) 
                break;
            ret_idx[*next_idx] = next_v_iter;
        }
    }
}

void compute_dip(const double *arr, const int *convex_m, int rel_length, int idx, int offset, double *ret_dip, int *ret_dip_idx) {
    assert(offset == 0 || offset == 1);
    
    int sign = (offset == 0) ? 1 : -1;
    *ret_dip = 0.;
    *ret_dip_idx = -1;  

    for (int j = idx; j < rel_length; ++j) {
        double temp_dip = 1.;
        int temp_dip_idx = -1;
        int  j_start = convex_m[j + 1 - offset], j_end = convex_m[j + offset];

        if (j_end - j_start > 1 && arr[j_end] != arr[j_start]) {

            double C = (j_end - j_start) / (arr[j_end] - arr[j_start]);

            for (int jj = j_start; jj <= j_end; ++jj) {
                double d = sign * (
                    (jj - j_start + sign) - (arr[jj] - arr[j_start]) * C
                );

                if (temp_dip < d) {
                    temp_dip = d; temp_dip_idx = jj;
                }
            }
        }

        if (*ret_dip < temp_dip) {
            *ret_dip = temp_dip; *ret_dip_idx = temp_dip_idx;
        }
        temp_dip = 1.; temp_dip_idx = -1;
    }
}

long double compute_largest_distance_greater_than_dip(const double *arr, const int *gcm, const int *lcm, int *ig, int *ix, int *ih, int *iv, int l_lcm, int debug) {
    /*	Find the largest distance greater than 'DIP' between the GCM and
        the LCM from LOW to HIGH.
    */
    long double ret_d = 0.;
    
    do {
        int gcm_ix = gcm[*ix], lcm_iv = lcm[*iv];
        int is_maj = gcm_ix > lcm_iv;
        int sign = (is_maj ? 1 : -1);

        int convex_m_ix = is_maj ? gcm_ix : lcm_iv;
        int convex_m_iv = is_maj ? lcm_iv : gcm_ix;
        int convex_m_i1 = is_maj ? gcm[*ix + 1] : lcm[*iv - 1];

        long double dx = sign * (
            (convex_m_iv - convex_m_i1 + sign) -
            ((long double) arr[convex_m_iv] - arr[convex_m_i1]) * 
                (convex_m_ix - convex_m_i1) / 
                (arr[convex_m_ix] - arr[convex_m_i1])
        );
        *iv = *iv + is_maj;
        *ix = *ix - (1 - is_maj);

        if (dx >= ret_d) {
		    ret_d = dx;
		    *ig = *ix + 1;
		    *ih = *iv - is_maj;
#if defined (DIPTEST_DEBUG)
		    if(debug >= 2) {
                cout << ((is_maj) ? "G" : "L") << "(" << (*ig) << ", " << (*ih) << ")";
            }
#endif // DIPTEST_DEBUG
	    }
        if (*ix < 1)	 *ix = 1;
	    if (*iv > l_lcm) *iv = l_lcm;

#if defined (DIPTEST_DEBUG)
	  if(debug) {
	      if (debug >= 2) {
            cout << " --> (ix, iv) = (" << (*ix) << ", " << (*iv) << ")" << endl;
        } else { 
            cout << ".";
        }
	  }
#endif // DIPTEST_DEBUG

    } while (gcm[*ix] != lcm[*iv]);

    return ret_d;
}

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
) {
/*
 * *low* contains the index of the current estimate  of the lower end.
 * of the modal interval, *high* contains the index for the upper end.
 * It can be used as: double xl = x[low], xu = x[high];
 */
#define low   lo_hi[0]
#define high  lo_hi[1]
/* 
 * *l_gcm* is defined as: relevant_length(GCM).
 * *l_lcm* is defined as: relevant_length(LCM)
 */
#define l_gcm lo_hi[2]
#define l_lcm lo_hi[3]

#ifndef DIPTEST_DEBUG
    UNUSED(debug);
#endif
    vector<int> range;
    double dip = (min_is_0) ? 0. : 1., dip_l, dip_u, tmp_max_dip;
    long double d = 0.; // TODO: check if with this 32-bit/64-bit differences go
    int ig, ih, iv, ix, i;
    int j_best, j_l, j_u;
    bool flag;

    /* Parameter adjustments, so that array referencing starts at 1,
       i.e., x[1]..x[n] 
    */
    --mj; --mn; --lcm; --gcm; --x;

    /* Consistency Checks */

    // non-positive check:
    if (n <= 0) {
        *ifault = 1;
        return 0.0;
    }
    
    // non-sorted array check:
    for (i = 2; i <= n; ++i)
        if (x[i] < x[i - 1]) {
            *ifault = 2;
            return 0.0;
        }

    low = 1;	high = n; 

    if (n < 2 || x[n] == x[1]) goto L_END;

#if defined (DIPTEST_DEBUG)
    if(debug)
	    cout << "'dip': START: (N = " << n << ")" 
             << " and 2N*dip = " << dip << "." 
             << endl;
#endif // DIPTEST_DEBUG

    /* Establish the indices   mn[1..n]  over which combination is necessary
       for the convex MINORANT (GCM) fit.
    */
    for( i = 1; i <= n; i++ )
        range.push_back(i);
    compute_convex_m_indices(x, mn, range);

    /* Establish the indices   mj[1..n]  over which combination is necessary
        for the concave MAJORANT (LCM) fit.
    */
    reverse(range.begin(), range.end());
    compute_convex_m_indices(x, mj, range);


    /* ------------------------- Start the cycling. ------------------------- */
    do {

        /* Collect the change points for the GCM from *high* to *low*. */
        gcm[1] = high;
        for(i = 1; gcm[i] > low; i++)
            gcm[i+1] = mn[gcm[i]];
        ig = l_gcm = i; // l_gcm == relevant_length(GCM)
        ix = ig - 1;   //  ix, ig  are counters for the convex minorant.

        /* Collect the change points for the LCM from *high* to *low*. */
        lcm[1] = low;
        for(i = 1; lcm[i] < high; i++)
            lcm[i+1] = mj[lcm[i]];
        ih = l_lcm = i; // l_lcm == relevant_length(LCM)
        iv = 2;        //  iv, ih  are counters for the concave majorant.

#if defined (DIPTEST_DEBUG)
        if(debug) {
            cout << "'dip': LOOP-BEGIN: 2n*D = " << dip
                 << " and [low, high] = ["
                 << setw(3) << low << ", " 
                 << setw(3) << high << "]"; 
            if(debug >= 3) {
                // Print the GCM:
                cout << " :" << endl << " gcm[1:" << l_gcm << "] = ";
                for(i = 1; i < l_gcm; i++)
                    cout << gcm[i] << ", ";
                cout << gcm[l_gcm] << endl;
                // Print the LCM:
                cout << " lcm[1:" << l_lcm << "] = ";
                for(i = 1; i < l_lcm; i++)
                    cout << lcm[i] << ", ";
                cout << lcm[l_lcm] << endl;
            } else {  // debug <= 2
                cout << "; (l_lcm, l_gcm) = (" 
                     << setw(2) << l_lcm << ", " 
                     << setw(2) << l_gcm << ")" << endl;
            }
        }
#endif // DIPTEST_DEBUG
        if (l_gcm != 2 || l_lcm != 2) {
#if defined (DIPTEST_DEBUG)
            if(debug) {
                cout << "'dip': CYCLE-BEGIN: while(gcm[ix] != lcm[iv])";
                if (debug >= 2) {
                    cout << endl;
                } else {
                    cout << " ";
                }
            }
#endif // DIPTEST_DEBUG

            d = compute_largest_distance_greater_than_dip(
                    x, gcm, lcm, &ig, &ix, &ih, &iv, l_lcm, debug
            );

#if defined (DIPTEST_DEBUG)      
            if(debug && debug < 2) cout << endl;
#endif // DIPTEST_DEBUG
        } else {

	        d = (min_is_0) ? 0. : 1.;

#if defined (DIPTEST_DEBUG)
	if(debug)
        cout << "'dip': NO-CYCLE: (l_lcm, l_gcm) = ("
             << setw(2) << l_lcm << ", " 
             << setw(2) << l_gcm << ") ==> d := "
             << d << endl;
#endif // DIPTEST_DEBUG
    }
    
    if (d < dip) goto L_END;

/*     Calculate the DIPs for the current LOW and HIGH. */
#if defined (DIPTEST_DEBUG)
    if(debug) cout << "'dip': MAIN-CALCULATION" << endl;
#endif // DIPTEST_DEBUG

    /* The DIP for the convex minorant. */
    compute_dip(x, gcm, l_gcm, ig, 0, &dip_l, &j_l);

    /* The DIP for the concave majorant. */
    compute_dip(x, lcm, l_lcm, ih, 1, &dip_u, &j_u);

#if defined (DIPTEST_DEBUG)
    if(debug)
        cout << " (dip_l, dip_u) = ("
             << dip_l << ", " 
             << dip_u << ")"; 
#endif // DIPTEST_DEBUG

    /* Determine the current maximum. */
    if(dip_u > dip_l) {
	    tmp_max_dip = dip_u; j_best = j_u;
    } else {
	    tmp_max_dip = dip_l; j_best = j_l;
    }
    if (dip < tmp_max_dip) {
	    dip = tmp_max_dip;
#if defined (DIPTEST_DEBUG)
	    if(debug)
            cout << " --> new larger dip " << tmp_max_dip
                 << " ( at index := " << j_best
                 << " )" << endl;
#endif // DIPTEST_DEBUG
    }

    flag = (low == gcm[ig] && high == lcm[ih]);
    low  = gcm[ig];
	high = lcm[ih];

#if defined (DIPTEST_DEBUG)
    if (flag && debug)
        cout << "'dip': LOOP-END: No improvement found neither in low := "
             << low << " nor in high := " << high << endl;
#endif // DIPTEST_DEBUG

} while (!flag);
/*---------------------------------------------------------------------------*/

L_END:
    /*  
     * M. Maechler -- speedup: Work with (2n * dip) everywhere but the very end!
     * It saves many divisions by n!
     */    
    dip /= (2*n);
    return dip;
} /* diptst */
#undef low
#undef high
#undef l_gcm
#undef l_lcm
