"""GWAS workflow adapters for shared RNA-seq inputs."""

from __future__ import annotations

from .rna_quantification import find_quantification_file

__all__ = ["find_quantification_file"]
