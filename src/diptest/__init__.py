from __future__ import annotations

import importlib.metadata

from diptest.diptest import dipstat, diptest
from diptest.lib import _diptest_core as _diptest
from diptest.lib._diptest_core import _has_openmp_support

__version__ = importlib.metadata.version("diptest")

__all__ = ["dipstat", "diptest", "_diptest", "_has_openmp_support", "__version__"]
