"""Amplicon-based metagenomics analysis (16S/ITS).

Provides OTU clustering, ASV denoising, and taxonomic classification
for marker gene sequencing data.
"""

from __future__ import annotations

from . import asv_denoising, otu_clustering, taxonomy

__all__ = [
    "otu_clustering",
    "asv_denoising",
    "taxonomy",
]
