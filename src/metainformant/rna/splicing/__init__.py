"""Alternative splicing analysis submodule for RNA-seq data.

This subpackage provides tools for detecting and quantifying alternative
splicing events from RNA-seq data, including:
- Splice junction detection and classification
- Percent Spliced In (PSI) computation
- Differential splicing analysis
- Novel junction identification
- Splice site strength scoring
- Isoform quantification via EM algorithm
- Splice graph construction and isoform enumeration
- Isoform diversity and usage comparison"""
from __future__ import annotations

from . import detection, isoforms, splice_analysis, splice_sites

__all__ = ['detection', 'isoforms', 'splice_analysis', 'splice_sites']
