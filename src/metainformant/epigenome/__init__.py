"""Epigenome domain functionality.

Public API:
- read_bedgraph: load bedGraph tracks into a DataFrame
- load_cpg_table / compute_beta_values / summarize_beta_by_chromosome: basic methylation utilities
"""

from .tracks import read_bedgraph  # noqa: F401
from .methylation import (  # noqa: F401
    load_cpg_table,
    compute_beta_values,
    summarize_beta_by_chromosome,
)

