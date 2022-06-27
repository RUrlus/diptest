import warnings
import psutil

import numpy as np

from diptest.lib import _diptest
from diptest.consts import Consts

_N_CORES = psutil.cpu_count(logical=False)
_N_CORES_MIN1 = _N_CORES - 1

if _diptest._has_openmp_support:
    _DEFAULT_N_THREADS = min(_N_CORES_MIN1, 4)
else:
    _DEFAULT_N_THREADS = 1


def dipstat(x, full_output=False, allow_zero=True, sort_x=True, debug=0):
    """
    Hartigan & Hartigan's dip statistic

    The dip statistic measures multimodality in a sample by the maximum
    difference, over all sample points, between the empirical distribution
    function, and the unimodal distribution function that minimizes that
    maximum difference.

    Parameters
    ----------
    x : np.ndarray
        the input samples
    full_output : boolean, default=False
        return dict alongside statistic, see below for details
    allow_zero : boolean, default=True
        if True the minimum value of the test statistic is
        allowed to be zero in cases where n <= 3 or all values in x
        are identical.
    sort_x : bool, default=True
        if False x is assumed to already be sorted in ascending order
    debug : int, default=0
        0 <= debug <= 3, print debugging messages, is ignored unless
        the pacakge was installed in debug mode

    Returns
    -------
    dip : double
        the dip statistic
    res : dict, optional
        returned if full_output == True.
        Contains the following fields:
            lo:     indices of lower end of modal interval
            hi:     indices of upper end of modal interval
            xl:     lower end of modal interval
            xu:     upper end of modal interval
            gcm:    (last-used) indices of the greatest concave majorant
            lcm:    (last-used) indices of the least concave majorant

    Reference
    -----------
    Hartigan, J. A., & Hartigan, P. M. (1985). The Dip Test of Unimodality.
        The Annals of Statistics.
    """
    if sort_x:
        x = np.sort(x)
    elif not isinstance(x, np.ndarray):
        x = np.asarray(x, order='C')
    elif not (x.flags.c_contiguous or x.flags.c_contiguous):
        x = np.copy(x, order='C')

    if ((x.ndim > 1) and not (
            (x.shape[1] == 1) or (x.shape[0]== 1)
        )
    ):
        raise TypeError("x should be one-dimensional")

    if full_output:
        res = _diptest.diptest_full(x, allow_zero, debug)
        dip = res.pop('dip')
        _gcm = res.pop('_gcm')
        res['gcm'] = _gcm[:res.pop('_lh_2')]
        _lcm = res.pop('_lcm')
        res['lcm'] = _lcm[:res.pop('_lh_3')]
        return float(dip), res
    return float(_diptest.diptest(x, allow_zero, debug))


def diptest(
    x,
    sort_x=True,
    allow_zero=True,
    boot_pval=False,
    n_boot=10000,
    n_threads=None,
    seed=None,
    stream=0,
):
    """
    Hartigan & Hartigan's dip test for unimodality.

    For X ~ F i.i.d., the null hypothesis is that F is a unimodal distribution.
    The alternative hypothesis is that F is multimodal (i.e. at least bimodal).
    Other than unimodality, the dip test does not assume any particular null
    distribution.

    Parameters
    ----------
    x : np.ndarray
        the input samples
    sort_x : bool, default=True
        if False x is assumed to already be sorted in ascending order
    allow_zero : boolean, default=True
        if True the minimum value of the test statistic is
        allowed to be zero in cases where n <= 3 or all values in x
        are identical.
    boot_pval : bool, default=False
        if True the p-value is computed using bootstrap samples from a
        uniform distribution, otherwise it is computed via linear
        interpolation of the tabulated critical values in dip_crit.txt.
    n_boot : int, default=10000
        if boot_pval=True, this sets the number of bootstrap samples to
        use for computing the p-value.
    n_threads : int, default=None
        number of threads to use when computing the p-value using bootstrap.
        Defaults to 4, if set to 1 the computation is
        performed single threaded. -1 will set the number of threads equal to
        all available cores
    seed : uint, default=None
        seed used for the generation of the uniform samples when computing the
        p-value.
    stream : uint, default=0
        stream used by PCG64_DXSM for a given seed. Note that setting the stream
        is only useful when you want to run multiple simulations based on a
        single seed. Parameter is ignored unless `boot_pval` is True and `seed`
        is not zero.

    Returns:
    -----------
    dip : double
        the dip statistic
    pval : double
        the p-value for the test

    Reference:
    -----------
    Hartigan, J. A., & Hartigan, P. M. (1985). The Dip Test of Unimodality.
        The Annals of Statistics.

    """
    n = x.size
    dip = dipstat(x, allow_zero=allow_zero, sort_x=sort_x)

    if n <= 3:
        warnings.warn('Dip test is not valid for n <= 3')
        return dip, 1.0

    if boot_pval:
        n_threads = n_threads or _DEFAULT_N_THREADS
        if n_threads == -1:
            n_threads = _N_CORES
        if n_threads > 1 and _diptest._has_openmp_support:
            pval = _diptest.diptest_pval_mt(
                dipstat=dip,
                n=n,
                n_boot=n_boot,
                allow_zero=allow_zero,
                seed=seed or 0,
                n_threads=n_threads
            )
            return dip, pval
        elif n_threads > 1 and not _diptest._has_openmp_support:
            warnings.warn("Extension was compiled without parallelisation support, ignoring ``n_threads``")
        if stream > Consts._UINT64_T_MAX:
            raise ValueError("`stream` must fit in a uint64_t.")
        pval = _diptest.diptest_pval(
            dipstat=dip,
            n=n,
            n_boot=n_boot,
            allow_zero=allow_zero,
            seed=seed or 0,
            stream=stream,
        )
        return dip, pval

    i1 = int(Consts._SAMPLE_SIZE.searchsorted(n, side='left'))
    i0 = i1 - 1

    # if n falls outside the range of tabulated sample sizes, use the
    # critical values for the nearest tabulated n (i.e. treat them as
    # 'asymptotic')
    i0 = max(0, i0)
    i1 = min(20, i1)

    # interpolate on sqrt(n)
    n0, n1 = Consts._SAMPLE_SIZE[[i0, i1]]

    y0 = np.sqrt(n0) * Consts._CRIT_VALS[i0]
    sD = np.sqrt(n) * dip
    if (i0 == i1):
        xp = y0
    else:
        fn = float(n - n0) / (n1 - n0)
        y1 = np.sqrt(n1) * Consts._CRIT_VALS[i1]
        xp = y0 + fn * (y1 - y0)

    pval = 1. - np.interp(sD, xp, Consts._ALPHA)

    return dip, float(pval)
