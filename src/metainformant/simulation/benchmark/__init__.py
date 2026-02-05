"""Benchmark dataset generation subpackage for simulation.

Provides synthetic benchmark dataset generators for classification,
regression, clustering, differential expression, GWAS, and method
comparison/evaluation.
"""

from __future__ import annotations

from .generators import (
    benchmark_suite,
    evaluate_benchmark,
    generate_benchmark_dataset,
    generate_synthetic_expression,
    generate_synthetic_variants,
)

__all__ = [
    "generate_benchmark_dataset",
    "generate_synthetic_variants",
    "generate_synthetic_expression",
    "evaluate_benchmark",
    "benchmark_suite",
]
