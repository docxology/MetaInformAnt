"""Shared utility functions for GWAS visualization.

Provides common helpers used across multiple visualization modules,
extracted to eliminate code duplication. Includes p-value key detection
and safe log10 transforms.
"""

from __future__ import annotations

import math
from typing import Any, Dict, Optional

# Floor value for p-value clamping to avoid log(0)
_P_VALUE_FLOOR: float = 1e-300


def detect_p_value_key(result_dict: Dict[str, Any]) -> Optional[str]:
    """Detect which key holds the p-value in a GWAS result dictionary.

    Checks for common p-value key names in priority order:
    p_value > pvalue > pval > P > p

    Args:
        result_dict: A single GWAS result dictionary.

    Returns:
        The key name holding the p-value, or None if no recognized key is found.
    """
    # Priority order matches the most common conventions across the codebase:
    # metainformant internal format, then PLINK/BOLT-LMM format, then minimal
    for key in ("p_value", "pvalue", "pval", "P", "p"):
        if key in result_dict:
            return key
    return None


def safe_log10_p(p_value: float) -> float:
    """Compute -log10(p) with clamping for zero and tiny values.

    Clamps p_value to [1e-300, +inf) before applying -log10, preventing
    -inf or math domain errors. Negative p-values are also clamped to
    the floor.

    Args:
        p_value: The p-value to transform. Should be in (0, 1] but
            zero, negative, and very small values are handled safely.

    Returns:
        -log10(clamped_p) as a float. Always finite and >= 0 for
        valid p-values.
    """
    clamped = max(p_value, _P_VALUE_FLOOR)
    return float(-math.log10(clamped))
