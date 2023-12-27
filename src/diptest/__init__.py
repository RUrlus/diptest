from __future__ import annotations

import importlib.metadata

from diptest.diptest import dipstat, diptest
from diptest.lib import _diptest_core as _diptest

__version__ = importlib.metadata.version("diptest")

__all__ = ["dipstat", "diptest", "_diptest", "__version__"]
