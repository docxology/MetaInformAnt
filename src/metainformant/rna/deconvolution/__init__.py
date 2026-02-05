"""Cell type deconvolution submodule for bulk RNA-seq data.

This subpackage provides tools for estimating cell type proportions
from bulk RNA-seq expression data, including:
- Non-negative least squares (NNLS) deconvolution
- SVR-based deconvolution (CIBERSORT-style)
- Signature matrix construction from reference profiles
- Marker gene selection
- Deconvolution result validation
- Batch deconvolution across multiple samples
"""

from __future__ import annotations

from .bulk_deconvolution import (
    batch_deconvolve,
    build_signature_matrix,
    deconvolve_nnls,
    deconvolve_svr,
    select_marker_genes,
    validate_deconvolution,
)

__all__ = [
    "deconvolve_nnls",
    "deconvolve_svr",
    "build_signature_matrix",
    "select_marker_genes",
    "validate_deconvolution",
    "batch_deconvolve",
]
