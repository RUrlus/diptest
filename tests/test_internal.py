import os
import warnings
import pytest
import numpy as np
import multiprocessing

import diptest as dt
from diptest.diptest import _mt_support

try:
    from diptest.lib._diptest import pcg_seed_test
    from diptest.lib._diptest import pcg_set_stream_test
    if _mt_support:
        from diptest.lib._diptest import pcg_mt_stream_test
    _has_support = True
except ImportError:
    _has_support = False

CORE_CNT = multiprocessing.cpu_count() - 1


if _has_support:

    def test_pcg_seed():
        N = 10000
        # seed == 0 triggers a drawn from random device and should result in
        # independent streams
        SEED = 0
        s0 = pcg_seed_test(N, SEED)
        s1 = pcg_seed_test(N, SEED)
        assert not np.isclose(s0, s1).all()

        SEED = 42
        s0 = pcg_seed_test(N, SEED)
        s1 = pcg_seed_test(N, SEED)
        assert np.isclose(s0, s1).all()

    def test_pcg_set_seed_stream():
        N = 100000
        streams = np.random.randint(
            0,
            np.iinfo(np.uint64).max,
            size=10,
            dtype=np.uint64
        )

        SEED = 42

        s0 = pcg_set_stream_test(size=N, seed=SEED, streams=streams)
        cov = np.cov(s0)
        cov[np.diag_indices(streams.size)] = 0
        assert (cov < 0.01).all()

        s1 = pcg_set_stream_test(size=N, seed=SEED, streams=streams)
        assert np.isclose(s0, s1).all()

        s1 = pcg_set_stream_test(size=N, seed=0, streams=streams)
        s2 = pcg_set_stream_test(size=N, seed=0, streams=streams)
        assert not np.isclose(s0, s1).all()
        assert not np.isclose(s1, s2).all()

    if _mt_support and CORE_CNT > 1:
        def test_pcg_mt_seed_stream():
            N = 100000
            SEED = 42

            s0 = pcg_mt_stream_test(row_size=N, seed=SEED, n_threads=CORE_CNT)
            cov = np.cov(s0)
            cov[np.diag_indices(CORE_CNT)] = 0
            assert (cov < 0.01).all()

            s1 = pcg_mt_stream_test(row_size=N, seed=SEED, n_threads=CORE_CNT)
            assert np.isclose(s0, s1).all()

            s1 = pcg_mt_stream_test(row_size=N, seed=0, n_threads=CORE_CNT)
            s2 = pcg_mt_stream_test(row_size=N, seed=0, n_threads=CORE_CNT)
            assert not np.isclose(s0, s1).all()
            assert not np.isclose(s1, s2).all()
