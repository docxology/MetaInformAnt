"""Life course and temporal phenotype analysis."""

from __future__ import annotations

from . import life_course, multivariate

try:
    from . import statistical
except ModuleNotFoundError as exc:  # pragma: no cover - depends on optional statsmodels availability
    if exc.name != "statsmodels":
        raise
    statistical = None  # type: ignore[assignment]

__all__ = ["life_course", "multivariate", "statistical"]
