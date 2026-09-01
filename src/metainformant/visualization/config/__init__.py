"""Configuration subpackage for visualization themes and palettes.

Re-exports all public APIs from :mod:`themes` and :mod:`palettes` for
convenience."""

from __future__ import annotations

import sys

from . import palettes, themes

# Self-reference so `from metainformant.visualization.config import config` works
config = sys.modules[__name__]

__all__ = ["conventions", "palettes", "themes", "config"]

from . import conventions  # noqa: E402
from .conventions import (  # noqa: E402,F401
    OKABE_ITO,
    add_log_zero_mark,
    apply_deterministic_rcparams,
    category_style,
    okabe_ito,
    save_figure_deterministic,
)
