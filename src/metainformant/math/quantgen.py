"""Quantitative genetics functions.

This module provides functions for quantitative genetics analysis,
including heritability estimation and selection response.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
import statistics

from metainformant.core import logging

logger = logging.get_logger(__name__)


def narrow_sense_heritability(additive_genetic_variance: float, phenotypic_variance: float) -> float:
    """Calculate narrow-sense heritability.

    Args:
        additive_genetic_variance: Additive genetic variance (VA)
        phenotypic_variance: Total phenotypic variance (VP)

    Returns:
        Narrow-sense heritability (h² = VA/VP)
    """
    if phenotypic_variance <= 0:
        raise ValueError("Phenotypic variance must be positive")
    if additive_genetic_variance < 0:
        raise ValueError("Additive genetic variance cannot be negative")

    return additive_genetic_variance / phenotypic_variance


def realized_heritability(selection_response: float, selection_differential: float) -> float:
    """Calculate realized heritability from selection experiment.

    Args:
        selection_response: Response to selection (R)
        selection_differential: Selection differential (S)

    Returns:
        Realized heritability (h²)
    """
    if selection_differential == 0:
        return 0.0

    return selection_response / selection_differential


def multilevel_selection_decomposition(group_trait_mean: float,
                                     individual_trait_variance: float,
                                     group_trait_variance: float,
                                     group_size: int) -> Dict[str, float]:
    """Decompose selection into individual and group-level components.

    Args:
        group_trait_mean: Mean trait value across groups
        individual_trait_variance: Variance within groups
        group_trait_variance: Variance between groups
        group_size: Average group size

    Returns:
        Dictionary with multilevel selection components
    """
    if group_size <= 0:
        raise ValueError("Group size must be positive")

    # Total phenotypic variance
    total_variance = individual_trait_variance + group_trait_variance * group_size

    # Individual-level selection
    individual_selection = individual_trait_variance / total_variance if total_variance > 0 else 0

    # Group-level selection
    group_selection = (group_trait_variance * group_size) / total_variance if total_variance > 0 else 0

    return {
        'total_variance': total_variance,
        'individual_level_selection': individual_selection,
        'group_level_selection': group_selection,
        'multilevel_selection_ratio': group_selection / individual_selection if individual_selection > 0 else 0
    }


def lande_equation_response(selection_intensity: float, heritability: float,
                          phenotypic_sd: float) -> float:
    """Calculate response to selection using Lande's equation.

    Args:
        selection_intensity: Selection intensity (i)
        heritability: Heritability (h²)
        phenotypic_sd: Phenotypic standard deviation

    Returns:
        Response to selection (R)
    """
    if phenotypic_sd <= 0:
        raise ValueError("Phenotypic standard deviation must be positive")

    # R = i * h² * σ_p
    return selection_intensity * heritability * phenotypic_sd


def heritability(broader_sense: bool, variance_components: Dict[str, float]) -> float:
    """Calculate heritability from variance components.

    Args:
        broader_sense: If True, calculate broad-sense heritability (H²), else narrow-sense (h²)
        variance_components: Dictionary with variance components

    Returns:
        Heritability estimate
    """
    required_keys = ['additive', 'dominance', 'epistasis', 'environmental']
    if not all(key in variance_components for key in required_keys):
        raise ValueError(f"Variance components must include: {required_keys}")

    additive = variance_components['additive']
    dominance = variance_components['dominance']
    epistasis = variance_components['epistasis']
    environmental = variance_components['environmental']

    total_genetic = additive + dominance + epistasis
    total_phenotypic = total_genetic + environmental

    if total_phenotypic == 0:
        return 0.0

    if broader_sense:
        # Broad-sense heritability: H² = Vg/Vp
        return total_genetic / total_phenotypic
    else:
        # Narrow-sense heritability: h² = Va/Vp
        return additive / total_phenotypic


def genetic_correlation(trait1_values: List[float], trait2_values: List[float],
                       relatedness_values: List[float]) -> float:
    """Calculate genetic correlation between traits.

    Args:
        trait1_values: Values of first trait
        trait2_values: Values of second trait
        relatedness_values: Relatedness coefficients

    Returns:
        Genetic correlation estimate
    """
    if len(trait1_values) != len(trait2_values) or len(trait1_values) != len(relatedness_values):
        raise ValueError("All input lists must have same length")

    if not trait1_values:
        return 0.0

    # Simplified estimation using cross-trait covariances
    # This is a basic approximation - full genetic correlation requires more complex estimation

    # Calculate trait means
    mean1 = statistics.mean(trait1_values)
    mean2 = statistics.mean(trait2_values)
    mean_r = statistics.mean(relatedness_values)

    # Calculate covariances weighted by relatedness
    cov_11 = sum(r * (t1 - mean1) ** 2 for t1, r in zip(trait1_values, relatedness_values)) / len(trait1_values)
    cov_22 = sum(r * (t2 - mean2) ** 2 for t2, r in zip(trait2_values, relatedness_values)) / len(trait2_values)
    cov_12 = sum(r * (t1 - mean1) * (t2 - mean2) for t1, t2, r in zip(trait1_values, trait2_values, relatedness_values)) / len(trait1_values)

    # Genetic correlation
    if cov_11 > 0 and cov_22 > 0:
        return cov_12 / (cov_11 ** 0.5 * cov_22 ** 0.5)
    else:
        return 0.0


def response_to_selection(selection_differential: float, heritability: float) -> float:
    """Calculate response to selection using breeder's equation.

    Args:
        selection_differential: Selection differential (S)
        heritability: Heritability (h²)

    Returns:
        Response to selection (R)
    """
    return heritability * selection_differential


def breeder_equation(selection_intensity: float, heritability: float,
                   phenotypic_sd: float) -> float:
    """Calculate response to selection using breeder's equation with selection intensity.

    Args:
        selection_intensity: Selection intensity (i)
        heritability: Heritability (h²)
        phenotypic_sd: Phenotypic standard deviation (σ_p)

    Returns:
        Response to selection (R = i * h² * σ_p)
    """
    if phenotypic_sd <= 0:
        raise ValueError("Phenotypic standard deviation must be positive")

    return selection_intensity * heritability * phenotypic_sd
