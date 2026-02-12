"""GWAS data retrieval and download utilities.

This module provides functions for downloading reference genomes, variant databases,
SRA sequencing data, and genome annotation for GWAS analysis."""
from __future__ import annotations

from . import config, download, genome, metadata, sra_download

__all__ = ['config', 'download', 'genome', 'metadata', 'sra_download']
