"""Population-genetics cross-domain analysis workflow.

Composes dna, gwas, math, and simulation domains into the population
genetics analysis pipeline; cross-domain composition is sanctioned here
(see ``scripts/quality/check_module_boundaries.py``).
"""

from __future__ import annotations

from .analysis import (
    analyze_dataset,
    compare_two_population_sequences,
    demographic_model_comparisons,
    genotype_structure_analysis,
    ld_summary,
    sequence_scenario_suite,
    summarize_scenario,
)
from .dataset_generation import generate_comprehensive_dataset
from .reporting import generate_summary_report
from .visualization import generate_visualizations

__all__ = [
    "analyze_dataset",
    "compare_two_population_sequences",
    "demographic_model_comparisons",
    "dataset_generation",
    "generate_comprehensive_dataset",
    "generate_summary_report",
    "generate_visualizations",
    "genotype_structure_analysis",
    "ld_summary",
    "reporting",
    "sequence_scenario_suite",
    "summarize_scenario",
    "visualization",
]
