"""Population genetics analysis module for METAINFORMANT.

Consolidates the population-genetics analysis surface: sequence-based
summary statistics, neutrality tests, two-population comparison, PCA /
kinship / HWE on genotype matrices, LD summaries, and demographic model
comparisons. The scenario analysis pipeline (as used by the synthetic
dataset workflow) is exposed via :func:`analyze_dataset`, with
``scripts/popgen/`` as its thin orchestrator.

Example:
    >>> from metainformant.popgen.workflow.analysis import summarize_scenario
    >>> result = summarize_scenario(
    ...     ["ATCG", "ATCG", "GCTA"], label="demo"
    ... )  # doctest: +SKIP
"""

from __future__ import annotations

from .workflow import analysis
from .workflow.analysis import (
    analyze_dataset,
    compare_two_population_sequences,
    demographic_model_comparisons,
    genotype_structure_analysis,
    ld_summary,
    sequence_scenario_suite,
    summarize_scenario,
)

__all__ = [
    "analysis",
    "analyze_dataset",
    "compare_two_population_sequences",
    "demographic_model_comparisons",
    "genotype_structure_analysis",
    "ld_summary",
    "sequence_scenario_suite",
    "summarize_scenario",
]
