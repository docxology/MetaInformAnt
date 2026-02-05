"""Psychophysical laws and contrast metrics.

This module implements fundamental laws of perception relating physical stimulus
intensity to perceived sensation.
"""

from __future__ import annotations

import math
from typing import Sequence, Union

import numpy as np


def weber_contrast(intensity: float, background: float) -> float:
    """Calculate Weber contrast.

    Suitable for small features on a uniform background.
    C_w = (I - I_b) / I_b

    Args:
        intensity: Stimulus intensity
        background: Background intensity

    Returns:
        Weber contrast value
    """
    if background == 0:
        return float("inf") if intensity != 0 else 0.0
    return (intensity - background) / background


def michelson_contrast(lt_max: float, lt_min: float) -> float:
    """Calculate Michelson contrast.

    Suitable for periodic patterns like gratings.
    C_m = (I_max - I_min) / (I_max + I_min)

    Args:
        lt_max: Maximum luminance
        lt_min: Minimum luminance

    Returns:
        Michelson contrast value (0.0 to 1.0)
    """
    if lt_max + lt_min == 0:
        return 0.0
    return (lt_max - lt_min) / (lt_max + lt_min)


def fechner_law(
    intensity: Union[float, np.ndarray], threshold: float = 1.0, k: float = 1.0
) -> Union[float, np.ndarray]:
    """Calculate sensation magnitude using Fechner's Law.

    S = k * ln(I / I_0)

    Args:
        intensity: Physical stimulus intensity (must be > 0)
        threshold: Absolute threshold of sensation (I_0)
        k: Weber fraction / Scaling constant

    Returns:
        Perceived sensation magnitude
    """
    intensity = np.asarray(intensity)
    # Avoid log(0) or log(negative)
    with np.errstate(invalid="ignore"):
        return k * np.log(intensity / threshold)


def stevens_power_law(intensity: Union[float, np.ndarray], exponent: float, k: float = 1.0) -> Union[float, np.ndarray]:
    """Calculate sensation magnitude using Stevens' Power Law.

    S = k * I^a

    Common exponents:
    - Brightness (point source): 0.5
    - Brightness (extended): 0.33
    - Loudness: 0.67
    - Electric shock: 3.5

    Args:
        intensity: Physical stimulus intensity
        exponent: Modality-specific exponent (a)
        k: Scaling constant

    Returns:
        Perceived sensation magnitude
    """
    return k * np.power(intensity, exponent)
