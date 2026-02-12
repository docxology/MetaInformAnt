"""GWAS integration subpackage for phenotype analysis.

Provides phenome-wide association studies (PheWAS), phenotype correlation,
genetic risk score computation, and heritability screening."""
from __future__ import annotations

from . import phewas

__all__ = ['phewas']
