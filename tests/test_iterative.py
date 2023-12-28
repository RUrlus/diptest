import numpy as np

import diptest as dt


def _ref_mirror_array(dat, left):
    """Mirrors dataset on input = [1, 2, 3], output = [-2, -1, 0, 1, 2]."""
    wdat = np.array(dat)
    if left:
        pivot = np.min(wdat)
        sdat = wdat - pivot
        mdat = np.concatenate((-sdat[sdat > 0], sdat))
    else:
        pivot = np.max(wdat)
        sdat = wdat - pivot
        mdat = np.concatenate((sdat, -sdat[sdat < 0]))
    return mdat


def test_mirror_array():
    assert np.array_equal(
        np.repeat(np.arange(-2, 3), 2).reshape(5, 2), dt.mirror_array(np.arange(1, 4))
    )
    assert dt.mirror_array(np.arange(10, dtype=np.int32)).dtype == np.int32
    assert dt.mirror_array(np.arange(10, dtype=np.int64)).dtype == np.int64
    assert dt.mirror_array(np.arange(10, dtype=np.float32)).dtype == np.float32
    assert dt.mirror_array(np.arange(10, dtype=np.float64)).dtype == np.float64

    rng = np.random.default_rng()
    sample = rng.normal(0, 1, 100)
    sort_sample = np.sort(sample)

    result = dt.mirror_array(sample)
    assert np.allclose(np.zeros(2), result[99, :])
    sort_result = dt.mirror_array(sort_sample)

    assert np.allclose(sort_result[:, 0], np.sort(_ref_mirror_array(sample, True)))
    assert np.allclose(sort_result[:, 1], np.sort(_ref_mirror_array(sample, False)))
