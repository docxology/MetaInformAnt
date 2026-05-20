"""Population genetics submodule."""

from __future__ import annotations

import math

from . import coalescent, core, demography, effective_size, fst, ld, selection, statistics
from .effective_size import effective_population_size_from_heterozygosity


def inbreeding_coefficient_from_fst(fst_value: float, subpopulations: int | None = None) -> float:
    """Estimate an inbreeding coefficient from F_ST."""
    if not (0 <= fst_value <= 1):
        raise ValueError("F_ST must be between 0 and 1")
    if fst_value == 1:
        return float("inf")
    return fst_value / (1.0 - fst_value)


def linkage_disequilibrium_decay_distance(r_squared: float, recombination_rate: float) -> float:
    """Distance at which LD decays to the supplied r² value."""
    if r_squared <= 0:
        return 0.0
    if r_squared >= 1 or recombination_rate <= 0:
        return float("inf")
    return -math.log(r_squared) / (2.0 * recombination_rate)


def coalescent_time_to_mrca(sample_size: int, effective_size: float) -> float:
    """Expected time to MRCA for a sample."""
    if sample_size <= 1:
        return 0.0
    return coalescent.expected_time_to_mrca(sample_size, effective_size)


__all__ = [
    "coalescent",
    "core",
    "demography",
    "effective_size",
    "fst",
    "ld",
    "selection",
    "statistics",
    "effective_population_size_from_heterozygosity",
    "inbreeding_coefficient_from_fst",
    "linkage_disequilibrium_decay_distance",
    "coalescent_time_to_mrca",
]
