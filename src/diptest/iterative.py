from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from diptest.lib import _diptest_core as _diptest


def mirror_array(
    arr: NDArray[np.int32 | int | np.float32 | float],
) -> NDArray[np.int32 | int | np.float32 | float]:
    """Mirror array both directions.

    The array is left mirrored (aka shift left/right s.t. min is zero and mirror around zero)
    and right mirrored (aka shift left/rigth s.t. the max is zero and mirror around zero).

    Args:
        arr: 1D array to be mirrored
    Returns:
        mirrors: 2D array of shape (2 * arr.size - 1, 2) where the first column contains the
        left mirror and the second column the right mirror

    """
    return _diptest.mirror_array(np.asarray(arr).ravel())
