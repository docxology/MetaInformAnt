"""Benchmark dataset generation subpackage for simulation.

Provides synthetic benchmark dataset generators for classification,
regression, clustering, differential expression, GWAS, and method
comparison/evaluation.
"""

from __future__ import annotations

from .generators import (
    generate_benchmark_dataset,
    generate_synthetic_variants,
    generate_synthetic_expression,
    evaluate_benchmark,
    benchmark_suite,
)

__all__ = [
    "generate_benchmark_dataset",
    "generate_synthetic_variants",
    "generate_synthetic_expression",
    "evaluate_benchmark",
    "benchmark_suite",
]
