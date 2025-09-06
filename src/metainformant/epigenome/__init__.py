"""Epigenome domain functionality.

Public API:
- read_bedgraph: load bedGraph tracks into a DataFrame
- load_cpg_table / compute_beta_values / summarize_beta_by_chromosome: basic methylation utilities
"""

from .methylation import compute_beta_values, load_cpg_table, summarize_beta_by_chromosome  # noqa: F401
from .tracks import read_bedgraph  # noqa: F401
