from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
import psutil

from diptest.consts import Consts
from diptest.lib import _diptest_core as _diptest

if TYPE_CHECKING:
    from typing import Any

_N_CORES = psutil.cpu_count(logical=False)
_N_CORES_MIN1 = _N_CORES - 1

_DEFAULT_N_THREADS = min(_N_CORES_MIN1, 4) if _diptest._has_openmp_support else 1


def dipstat(
    x, full_output=False, allow_zero=True, sort_x=True, debug=0
) -> float | tuple[float, dict[str, Any]]:
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
        the package was installed in debug mode

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
            dipidx: index of the dip

    Reference
    -----------
    Hartigan, J. A., & Hartigan, P. M. (1985). The Dip Test of Unimodality.
        The Annals of Statistics.
    """
    if sort_x:
        x = np.sort(x)
    elif not isinstance(x, np.ndarray):
        x = np.asarray(x, order="C")
    elif not (x.flags.c_contiguous or x.flags.f_contiguous):
        x = np.copy(x, order="C")

    if (x.ndim > 1) and not ((x.shape[1] == 1) or (x.shape[0] == 1)):
        msg = "x should be one-dimensional"
        raise TypeError(msg)

    if full_output:
        res = _diptest.diptest_full(x, allow_zero, debug)
        dip = res.pop("dip")
        _gcm = res.pop("_gcm")
        res["gcm"] = _gcm[: res.pop("_lh_2")]
        _lcm = res.pop("_lcm")
        res["lcm"] = _lcm[: res.pop("_lh_3")]
        return float(dip), res
    return float(_diptest.diptest(x, allow_zero, debug))


def diptest(
    x,
    full_output=False,
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
    full_output : boolean, default=False
        return dict alongside statistic, see below for details
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
    res : dict, optional
        returned if full_output == True.
        Contains the following fields:
            lo:     indices of lower end of modal interval
            hi:     indices of upper end of modal interval
            xl:     lower end of modal interval
            xu:     upper end of modal interval
            gcm:    (last-used) indices of the greatest concave majorant
            lcm:    (last-used) indices of the least concave majorant
            dipidx: index of the dip

    Reference:
    -----------
    Hartigan, J. A., & Hartigan, P. M. (1985). The Dip Test of Unimodality.
        The Annals of Statistics.

    """
    n = x.size
    r = dipstat(x, full_output=full_output, allow_zero=allow_zero, sort_x=sort_x)
    dip = r if not full_output else r[0]

    if n <= 3:
        warnings.warn("Dip test is not valid for n <= 3", stacklevel=1)
        return (dip, 1.0) if not full_output else (dip, 1.0, r[1])

    if boot_pval:
        n_threads = n_threads or _DEFAULT_N_THREADS
        if n_threads == -1:
            n_threads = _N_CORES

        kwargs = {
            "dipstat": dip,
            "n": n,
            "n_boot": n_boot,
            "allow_zero": allow_zero,
            "seed": seed or 0,
        }

        func = None
        if n_threads > 1:
            if _diptest._has_openmp_support:
                kwargs["n_threads"] = n_threads
                func = _diptest.diptest_pval_mt
            else:
                msg = "Extension was compiled without parallelisation support, ignoring ``n_threads``"
                warnings.warn(msg, stacklevel=1)

        if func is None:
            if stream > Consts._UINT64_T_MAX:
                msg = "`stream` must fit in a uint64_t."
                raise ValueError(msg)
            kwargs["stream"] = stream
            func = _diptest.diptest_pval

        pval = func(**kwargs)
    else:
        pval = Consts.compute_pval_interpolation(n, dip)

    if full_output:
        return dip, pval, r[1]
    return dip, pval
