"""Cell type deconvolution submodule for bulk RNA-seq data.

This subpackage provides tools for estimating cell type proportions
from bulk RNA-seq expression data, including:
- Non-negative least squares (NNLS) deconvolution
- SVR-based deconvolution (CIBERSORT-style)
- Signature matrix construction from reference profiles
- Marker gene selection
- Deconvolution result validation
- Batch deconvolution across multiple samples"""
from __future__ import annotations

from . import bulk_deconvolution

__all__ = ['bulk_deconvolution']
