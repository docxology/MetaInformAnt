"""GWAS visualization modules.

Organized into subpackages:
- population/: Population structure and geographic visualizations
- statistical/: Statistical analysis and comparison plots
- genomic/: Genome-wide, regional, LD, and variant visualizations
- interactive/: Interactive plots, composite panels, fine-mapping, phenotype, suite"""
from __future__ import annotations

from . import config, eqtl_visualization, general, genomic, interactive, population, statistical, utils
from .general import manhattan_plot, qq_plot
from .utils import detect_p_value_key, safe_log10_p

__all__ = [
    'config', 'eqtl_visualization', 'general', 'genomic', 'interactive',
    'population', 'statistical', 'utils',
    'manhattan_plot', 'qq_plot', 'detect_p_value_key', 'safe_log10_p',
]
