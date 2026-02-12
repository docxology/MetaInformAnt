"""Protein function prediction utilities.

This module provides tools for predicting protein function from domain
composition, subcellular localization, solubility, physicochemical properties,
intrinsically disordered regions, active sites, and post-translational
modification sites. All predictions use pure Python implementations with
biologically grounded algorithms and scoring matrices.

This module re-exports all public symbols from :mod:`prediction_core` and
:mod:`prediction_analysis` for backward compatibility.
"""

from __future__ import annotations

from metainformant.protein.function.prediction_core import (
    _AA_MW,
    _DISORDER_PROPENSITY,
    _DOMAIN_FUNCTION_MAP,
    _ER_RETENTION_PATTERNS,
    _KD_HYDROPHOBICITY,
    _MTS_PATTERN,
    _NLS_PATTERNS,
    _PEROXISOMAL_SIGNALS,
    _PKA_CTERM,
    _PKA_NTERM,
    _PKA_SIDE,
    _charge_at_ph,
    _compute_composition,
    _compute_instability_index,
    _compute_isoelectric_point,
    compute_physicochemical,
    predict_function_from_domains,
    predict_localization,
    predict_solubility,
)
from metainformant.protein.function.prediction_analysis import (
    _ACTIVE_SITE_PATTERNS,
    _PTM_PATTERNS,
    find_active_sites,
    predict_disordered_regions,
    predict_post_translational_mods,
)

__all__ = [
    # Core prediction functions
    "predict_function_from_domains",
    "predict_localization",
    "predict_solubility",
    "compute_physicochemical",
    # Analysis functions
    "predict_disordered_regions",
    "find_active_sites",
    "predict_post_translational_mods",
    # Helper functions (private but re-exported for compatibility)
    "_compute_composition",
    "_compute_isoelectric_point",
    "_charge_at_ph",
    "_compute_instability_index",
    # Constants (private but re-exported for compatibility)
    "_KD_HYDROPHOBICITY",
    "_AA_MW",
    "_PKA_NTERM",
    "_PKA_CTERM",
    "_PKA_SIDE",
    "_DISORDER_PROPENSITY",
    "_DOMAIN_FUNCTION_MAP",
    "_NLS_PATTERNS",
    "_MTS_PATTERN",
    "_ER_RETENTION_PATTERNS",
    "_PEROXISOMAL_SIGNALS",
    "_ACTIVE_SITE_PATTERNS",
    "_PTM_PATTERNS",
]
