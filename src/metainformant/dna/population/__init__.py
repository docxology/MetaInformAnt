"""Population genetics analysis sub-package."""
from __future__ import annotations

from . import analysis, core, visualization, visualization_core, visualization_stats

# Re-export public functions from core.py
from .core import (
    allele_frequencies,
    expected_heterozygosity,
    fay_wu_h_from_sequences,
    fixation_index,
    fu_and_li_d_star_from_sequences,
    fu_and_li_f_star_from_sequences,
    hardy_weinberg_allele_freqs,
    hudson_fst,
    linkage_disequilibrium,
    nucleotide_diversity,
    observed_heterozygosity,
    segregating_sites,
    tajimas_d,
    wattersons_theta,
)

# Re-export public functions from analysis.py
from .analysis import (
    calculate_fst,
    calculate_fu_li_d,
    calculate_ld_decay,
    calculate_nucleotide_diversity,
    calculate_segregating_sites,
    calculate_summary_statistics,
    calculate_tajima_d,
    calculate_wattersons_theta,
    detect_population_structure,
    detect_selection,
    estimate_population_size,
    mcdonald_kreitman_test,
)

__all__ = [
    "analysis", "core", "visualization", "visualization_core", "visualization_stats",
    # core.py
    "allele_frequencies", "expected_heterozygosity", "fay_wu_h_from_sequences",
    "fixation_index", "fu_and_li_d_star_from_sequences", "fu_and_li_f_star_from_sequences",
    "hardy_weinberg_allele_freqs", "hudson_fst", "linkage_disequilibrium",
    "nucleotide_diversity", "observed_heterozygosity", "segregating_sites",
    "tajimas_d", "wattersons_theta",
    # analysis.py
    "calculate_fst", "calculate_fu_li_d", "calculate_ld_decay",
    "calculate_nucleotide_diversity", "calculate_segregating_sites",
    "calculate_summary_statistics", "calculate_tajima_d", "calculate_wattersons_theta",
    "detect_population_structure", "detect_selection", "estimate_population_size",
    "mcdonald_kreitman_test",
]

