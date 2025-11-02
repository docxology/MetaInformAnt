from __future__ import annotations

from typing import Sequence, Tuple


def narrow_sense_heritability(additive_variance: float, phenotypic_variance: float) -> float:
    """Calculate narrow-sense heritability h² = VA / VP.
    
    Narrow-sense heritability is the proportion of phenotypic variance
    explained by additive genetic variance. This is the component that
    responds to selection.
    
    Args:
        additive_variance: Additive genetic variance (VA)
        phenotypic_variance: Total phenotypic variance (VP)
        
    Returns:
        Narrow-sense heritability in [0, 1]. Returns 0.0 if VP <= 0.
        Values clamped to [0, 1] range.
        
    Examples:
        >>> narrow_sense_heritability(additive_variance=0.3, phenotypic_variance=0.5)
        0.6
        >>> narrow_sense_heritability(additive_variance=0.2, phenotypic_variance=1.0)
        0.2
        
    References:
        Falconer, D. S., & Mackay, T. F. C. (1996). Introduction to quantitative
        genetics. Longman.
    """
    VA = max(0.0, additive_variance)
    VP = max(0.0, phenotypic_variance)
    if VP <= 0:
        return 0.0
    h2 = VA / VP
    if h2 < 0:
        return 0.0
    if h2 > 1:
        return 1.0
    return h2


def breeders_equation_response(selection_differential: float, heritability: float) -> float:
    """Calculate response to selection using breeder's equation: R = h² × S.
    
    The breeder's equation predicts the change in mean phenotype after
    one generation of selection, given selection differential and heritability.
    
    Args:
        selection_differential: Difference between mean of selected parents
            and population mean (S)
        heritability: Narrow-sense heritability (h²), clamped to [0, 1]
            
    Returns:
        Expected response to selection (R)
        
    Examples:
        >>> breeders_equation_response(selection_differential=2.0, heritability=0.5)
        1.0
        >>> breeders_equation_response(selection_differential=1.5, heritability=0.8)
        1.2
        
    References:
        Lush, J. L. (1937). Animal breeding plans. Iowa State College Press.
    """
    h2 = max(0.0, min(1.0, heritability))
    return h2 * selection_differential


def lande_equation_response(gradient: Sequence[float], G_matrix: Sequence[Sequence[float]]) -> Tuple[float, ...]:
    """Calculate multivariate response to selection: Δz = G β.
    
    Lande's equation extends the breeder's equation to multiple traits,
    accounting for genetic correlations between traits through the G matrix.
    
    Args:
        gradient: Selection gradient vector (β), one element per trait
        G_matrix: Genetic variance-covariance matrix (G). Should be symmetric
            positive semidefinite. Each row/column corresponds to a trait.
            
    Returns:
        Tuple of response vectors (Δz), one element per trait.
        Returns empty tuple if dimensions don't match or inputs invalid.
        
    Examples:
        >>> gradient = [0.5, 0.3]
        >>> G = [[0.4, 0.1], [0.1, 0.3]]
        >>> response = lande_equation_response(gradient, G)
        >>> len(response)
        2
        
    References:
        Lande, R. (1979). Quantitative genetic analysis of multivariate
        evolution, applied to brain:body size allometry. Evolution, 33(1), 402-416.
    """
    beta = list(gradient)
    G = [list(row) for row in G_matrix]
    if not G or not beta or any(len(row) != len(beta) for row in G):
        return tuple()
    response = [sum(G[i][j] * beta[j] for j in range(len(beta))) for i in range(len(beta))]
    return tuple(response)


def realized_heritability(response: float, selection_differential: float) -> float:
    """h^2_realized = R / S (guarding S=0)."""
    if abs(selection_differential) < 1e-12:
        return 0.0
    h2 = response / selection_differential
    return max(0.0, min(1.0, h2))
