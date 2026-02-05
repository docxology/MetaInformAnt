"""Community ecology analysis and biodiversity metrics.

This module provides comprehensive tools for analyzing ecological communities,
including diversity indices, species abundance distributions, community
similarity measures, and biodiversity assessment.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple

from metainformant.core import errors, logging, validation

logger = logging.get_logger(__name__)


def calculate_diversity(
    species_matrix: List[List[float]] | Dict[str, List[float]], method: str = "shannon"
) -> List[float] | Dict[str, List[float]]:
    """Calculate diversity indices for ecological communities.

    Args:
        species_matrix: Species abundance data. Can be:
            - List of lists: Each inner list represents a community
            - Dict: Keys are community names, values are abundance lists
        method: Diversity index to calculate ("shannon", "simpson", "invsimpson", "richness")

    Returns:
        Diversity values for each community (same structure as input)

    Raises:
        ValueError: If method is not supported or data format is invalid
    """
    validation.validate_not_empty(species_matrix, "species_matrix")
    supported_methods = ["shannon", "simpson", "invsimpson", "richness"]
    if method not in supported_methods:
        raise ValueError(f"Unsupported diversity method: {method}. Choose from {supported_methods}")

    logger.info(f"Calculating {method} diversity for ecological communities")

    if isinstance(species_matrix, dict):
        results = {}
        for community_name, abundances in species_matrix.items():
            results[community_name] = [calculate_single_diversity(abundances, method)]
        return results
    else:
        results = []
        for abundances in species_matrix:
            results.append(calculate_single_diversity(abundances, method))
        return results


def calculate_single_diversity(abundances: List[float], method: str) -> float:
    """Calculate diversity index for a single community.

    Args:
        abundances: List of species abundances
        method: Diversity calculation method

    Returns:
        Diversity index value
    """
    # Filter out zero abundances and ensure positive values
    abundances = [max(0, x) for x in abundances if x > 0]

    if not abundances:
        return 0.0

    total_abundance = sum(abundances)

    if total_abundance == 0:
        return 0.0

    # Convert to proportions
    proportions = [abundance / total_abundance for abundance in abundances]

    if method == "shannon":
        # Shannon diversity index: H = -Σ(pi * ln(pi))
        return -sum(p * math.log(p) for p in proportions if p > 0)

    elif method == "simpson":
        # Simpson diversity index: 1 - D, where D = Σ(pi^2) is the probability of picking two individuals of the same species
        simpson_concentration = sum(p**2 for p in proportions)
        return 1.0 - simpson_concentration

    elif method == "invsimpson":
        # Inverse Simpson diversity index: 1/D
        simpson_d = sum(p**2 for p in proportions)
        return 1.0 / simpson_d if simpson_d > 0 else 0.0

    elif method == "richness":
        # Species richness: number of species present
        return len(abundances)

    else:
        raise ValueError(f"Unknown diversity method: {method}")


def species_richness(
    community_data: List[float] | List[List[float]] | Dict[str, List[float]],
) -> int | List[int] | Dict[str, int]:
    """Calculate species richness for one or more communities.

    Species richness is the number of different species present.

    Args:
        community_data: Species abundance data. Can be:
            - Simple list of floats: returns int count
            - List of lists: returns list of int counts
            - Dict of lists: returns dict of int counts

    Returns:
        Species richness count(s)
    """
    if not community_data:
        return 0

    # Simple list case (flat list of numbers, not nested)
    if isinstance(community_data, list) and (not community_data or not isinstance(community_data[0], list)):
        return species_richness_simple(community_data)

    # Dict case
    if isinstance(community_data, dict):
        results = {}
        for community_name, abundances in community_data.items():
            richness = sum(1 for abundance in abundances if abundance > 0)
            results[community_name] = richness
        return results

    # Nested list case
    results = []
    for abundances in community_data:
        richness = sum(1 for abundance in abundances if abundance > 0)
        results.append(richness)
    return results


def calculate_evenness(abundances: List[float], method: str = "pielou") -> float:
    """Calculate species evenness (equitability) for a community.

    Args:
        abundances: Species abundance data
        method: Evenness calculation method ("pielou", "simpson")

    Returns:
        Evenness index value (0-1, where 1 is perfect evenness)

    Raises:
        ValueError: If method is not supported
    """
    if not abundances:
        return 0.0

    supported_methods = ["pielou", "simpson"]
    if method not in supported_methods:
        raise ValueError(f"Unsupported evenness method: {method}. Choose from {supported_methods}")

    # Calculate Shannon diversity
    shannon_h = calculate_single_diversity(abundances, "shannon")

    # Filter species with positive abundance
    positive_abundances = [x for x in abundances if x > 0]
    richness = len(positive_abundances)

    if richness <= 1:
        return 0.0  # Cannot have evenness with 0 or 1 species

    if method == "pielou":
        # Pielou's evenness: J = H / ln(S)
        return shannon_h / math.log(richness)

    elif method == "simpson":
        # Simpson's evenness: E = 1/D / S = (1/D) * (1/S)
        simpson_d = calculate_single_diversity(abundances, "simpson")
        if simpson_d > 0:
            return (1.0 / simpson_d) / richness
        else:
            return 0.0

    return 0.0


def rarefaction_curve(abundances: List[float], max_samples: Optional[int] = None) -> List[Tuple[int, float]]:
    """Calculate rarefaction curve for species richness estimation.

    Args:
        abundances: Species abundance data
        max_samples: Maximum number of samples to evaluate

    Returns:
        List of (sample_size, richness) tuples
    """
    # Filter and sort abundances (descending)
    abundances = sorted([int(x) for x in abundances if x > 0], reverse=True)

    if not abundances:
        return [(0, 0.0)]

    total_individuals = sum(abundances)

    if max_samples is None:
        max_samples = min(total_individuals, 1000)  # Default max

    max_samples = min(max_samples, total_individuals)

    curve = []

    for n in range(1, max_samples + 1):
        # Calculate expected richness for sample size n
        expected_richness = 0.0

        for abundance in abundances:
            if abundance > 0:
                # Hypergeometric expectation
                prob = 1.0 - math.exp(-n * abundance / total_individuals)
                expected_richness += prob

        curve.append((n, expected_richness))

    return curve


def species_accumulation_curve(sampling_effort: List[int]) -> List[Tuple[int, float]]:
    """Calculate species accumulation curve from sampling effort data.

    Args:
        sampling_effort: List of cumulative species counts at each sampling effort

    Returns:
        Smoothed accumulation curve as (effort, richness) tuples
    """
    if not sampling_effort:
        return [(0, 0.0)]

    curve = []
    for i, richness in enumerate(sampling_effort):
        effort = i + 1
        curve.append((effort, float(richness)))

    return curve


def beta_diversity(community1: List[float], community2: List[float], method: str = "bray_curtis") -> float:
    """Calculate beta diversity between two communities.

    Beta diversity measures the difference in species composition between communities.

    Args:
        community1: Species abundances for first community
        community2: Species abundances for second community
        method: Beta diversity method ("bray_curtis", "jaccard", "sorensen")

    Returns:
        Beta diversity index (0 = identical, higher = more different)

    Raises:
        ValueError: If method is not supported
    """
    supported_methods = ["bray_curtis", "jaccard", "sorensen"]
    if method not in supported_methods:
        raise ValueError(f"Unsupported beta diversity method: {method}. Choose from {supported_methods}")

    # Ensure same length (pad with zeros if necessary)
    max_len = max(len(community1), len(community2))
    comm1 = community1 + [0] * (max_len - len(community1))
    comm2 = community2 + [0] * (max_len - len(community2))

    if method == "bray_curtis":
        # Bray-Curtis dissimilarity: sum|xi - yi| / sum(xi + yi)
        numerator = sum(abs(a - b) for a, b in zip(comm1, comm2))
        denominator = sum(a + b for a, b in zip(comm1, comm2))

        return numerator / denominator if denominator > 0 else 0.0

    elif method == "jaccard":
        # Jaccard distance: 1 - (intersection / union)
        set1 = set(i for i, v in enumerate(comm1) if v > 0)
        set2 = set(i for i, v in enumerate(comm2) if v > 0)

        intersection = len(set1.intersection(set2))
        union = len(set1.union(set2))

        if union == 0:
            return 0.0

        return 1.0 - (intersection / union)

    elif method == "sorensen":
        # Sørensen distance: 1 - (2 * intersection / (sum1 + sum2))
        set1 = set(i for i, v in enumerate(comm1) if v > 0)
        set2 = set(i for i, v in enumerate(comm2) if v > 0)

        intersection = len(set1.intersection(set2))
        sum_sets = len(set1) + len(set2)

        if sum_sets == 0:
            return 0.0

        return 1.0 - (2 * intersection / sum_sets)

    return 0.0


def rank_abundance_curve(abundances: List[float]) -> List[Tuple[int, float]]:
    """Generate rank-abundance curve (Whittaker plot).

    Args:
        abundances: Species abundance data

    Returns:
        List of (rank, abundance) tuples sorted by abundance (descending)
    """
    # Sort abundances in descending order
    sorted_abundances = sorted(abundances, reverse=True)

    curve = []
    for rank, abundance in enumerate(sorted_abundances, 1):
        curve.append((rank, abundance))

    return curve


def dominance_diversity_curve(abundances: List[float]) -> List[Tuple[float, float]]:
    """Calculate dominance-diversity curve (k-dominance curve).

    Args:
        abundances: Species abundance data

    Returns:
        List of (dominance_level, diversity) tuples
    """
    if not abundances:
        return [(0.0, 0.0)]

    # Sort in descending order
    sorted_abundances = sorted(abundances, reverse=True)
    total_abundance = sum(sorted_abundances)

    if total_abundance == 0:
        return [(0.0, 0.0)]

    curve = []
    cumulative_abundance = 0.0

    for i, abundance in enumerate(sorted_abundances):
        cumulative_abundance += abundance
        dominance = cumulative_abundance / total_abundance

        # Number of species needed to reach this dominance level
        diversity = i + 1

        curve.append((dominance, diversity))

    return curve


def species_area_relationship(species_counts: List[int], area_sizes: List[float]) -> Dict[str, Any]:
    """Fit species-area relationship (Arrhenius/power law model).

    Args:
        species_counts: Number of species at each area
        area_sizes: Corresponding area sizes

    Returns:
        Dictionary with fitted parameters and statistics
    """
    if len(species_counts) != len(area_sizes) or len(species_counts) < 2:
        raise ValueError("Need at least 2 corresponding species counts and area sizes")

    # Log-transform data for linear regression
    log_areas = [math.log(area) for area in area_sizes]
    log_species = [math.log(count) for count in species_counts]

    # Linear regression: log(S) = log(c) + z * log(A)
    try:
        result = statistics.linear_regression(log_areas, log_species)
        slope = result.slope
        intercept = result.intercept

        # Convert back to power law parameters: S = c * A^z
        c = math.exp(intercept)
        z = slope

        # Calculate R-squared manually
        mean_y = statistics.mean(log_species)
        ss_tot = sum((y - mean_y) ** 2 for y in log_species)
        ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(log_areas, log_species))
        r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

        return {
            "c": c,
            "z": z,
            "r_squared": r_squared,
            "intercept": intercept,
            "slope": slope,
            "model": "power_law",
            "equation": f"S = {c:.4f} * A^{z:.4f}",
        }

    except Exception as e:
        logger.warning(f"Could not fit species-area relationship: {e}")
        return {
            "error": str(e),
            "c": None,
            "z": None,
            "r_squared": None,
        }


def nestedness_temperature_calculator(presence_absence_matrix: List[List[int]]) -> float:
    """Calculate nestedness temperature (NTC) for a community matrix.

    Nestedness measures how much the species composition of smaller communities
    is a subset of larger communities.

    Args:
        presence_absence_matrix: Binary matrix (rows = sites, columns = species)

    Returns:
        Nestedness temperature (0 = perfectly nested, 100 = random)
    """
    if not presence_absence_matrix:
        return 100.0

    # Convert to numpy-like operations (simplified implementation)
    matrix = presence_absence_matrix
    n_sites = len(matrix)
    n_species = len(matrix[0]) if matrix else 0

    if n_sites < 2 or n_species < 2:
        return 100.0

    # Count unexpected presences and absences
    unexpected = 0
    total_comparisons = 0

    # Sort sites by richness (descending)
    site_richness = [(i, sum(row)) for i, row in enumerate(matrix)]
    site_richness.sort(key=lambda x: x[1], reverse=True)

    for i in range(n_sites):
        for j in range(i + 1, n_sites):
            site_i = site_richness[i][0]
            site_j = site_richness[j][0]

            # Compare species presence
            for species in range(n_species):
                present_i = matrix[site_i][species]
                present_j = matrix[site_j][species]

                # If species is in poorer site but not in richer site = unexpected
                if present_j and not present_i:
                    unexpected += 1

                total_comparisons += 1

    # Calculate temperature
    if total_comparisons == 0:
        return 100.0

    temperature = (unexpected / total_comparisons) * 100
    return temperature


def calculate_biodiversity_indices(
    community_data: List[List[float]], indices: Optional[List[str]] = None
) -> Dict[str, List[float]]:
    """Calculate multiple biodiversity indices for communities.

    Args:
        community_data: List of community abundance vectors
        indices: List of indices to calculate (default: all available)

    Returns:
        Dictionary mapping index names to lists of values
    """
    if indices is None:
        indices = ["shannon", "simpson", "invsimpson", "richness"]

    results = {}

    for index_name in indices:
        try:
            values = calculate_diversity(community_data, index_name)
            if isinstance(values, list):
                results[index_name] = values
            else:
                results[index_name] = [values]
        except Exception as e:
            logger.warning(f"Could not calculate {index_name} index: {e}")
            results[index_name] = [0.0] * len(community_data)

    return results


def community_similarity_matrix(communities: List[List[float]], method: str = "bray_curtis") -> List[List[float]]:
    """Calculate pairwise similarity matrix for communities.

    Args:
        communities: List of community abundance vectors
        method: Similarity/dissimilarity method

    Returns:
        Symmetric similarity matrix
    """
    n_communities = len(communities)
    similarity_matrix = [[0.0] * n_communities for _ in range(n_communities)]

    for i in range(n_communities):
        for j in range(i, n_communities):
            if method in ["bray_curtis", "jaccard", "sorensen"]:
                # These are dissimilarity measures, convert to similarity
                dissimilarity = beta_diversity(communities[i], communities[j], method)
                similarity = 1.0 - dissimilarity
            else:
                # Assume similarity measure
                similarity = beta_diversity(communities[i], communities[j], method)

            similarity_matrix[i][j] = similarity
            similarity_matrix[j][i] = similarity

    return similarity_matrix


def alpha_beta_gamma_diversity(communities: List[List[float]]) -> Dict[str, float]:
    """Calculate alpha, beta, and gamma diversity for a set of communities.

    Alpha: Average diversity within communities
    Beta: Diversity between communities
    Gamma: Total diversity across all communities

    Args:
        communities: List of community abundance vectors

    Returns:
        Dictionary with alpha, beta, and gamma diversity values
    """
    if not communities:
        return {"alpha": 0.0, "beta": 0.0, "gamma": 0.0}

    # Alpha diversity (mean Shannon diversity)
    alpha_values = calculate_diversity(communities, "shannon")
    if isinstance(alpha_values, list):
        alpha = statistics.mean(alpha_values) if alpha_values else 0.0
    else:
        alpha = alpha_values

    # Gamma diversity (diversity of pooled communities)
    pooled_abundances = []
    max_species = max(len(comm) for comm in communities) if communities else 0

    for species_idx in range(max_species):
        total_abundance = sum(comm[species_idx] if species_idx < len(comm) else 0 for comm in communities)
        pooled_abundances.append(total_abundance)

    gamma = calculate_single_diversity(pooled_abundances, "shannon")

    # Beta diversity (gamma / alpha)
    beta = gamma / alpha if alpha > 0 else 0.0

    return {
        "alpha": alpha,
        "beta": beta,
        "gamma": gamma,
    }


def generate_ecology_report(
    community_data: List[List[float]],
    sample_names: Optional[List[str]] = None,
    output_path: Optional[str | Path] = None,
) -> str:
    """Generate a comprehensive ecology analysis report.

    Args:
        community_data: List of community abundance vectors
        sample_names: Optional names for communities
        output_path: Optional path to save the report

    Returns:
        Formatted ecology report
    """
    if sample_names is None:
        sample_names = [f"Community_{i+1}" for i in range(len(community_data))]

    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("ECOLOGICAL COMMUNITY ANALYSIS REPORT")
    report_lines.append("=" * 60)
    report_lines.append("")

    # Basic statistics
    report_lines.append(f"Number of Communities: {len(community_data)}")
    total_species = len(set(sum([list(range(len(comm))) for comm in community_data], [])))
    report_lines.append(f"Total Species Detected: {total_species}")
    report_lines.append("")

    # Diversity indices
    diversity_indices = calculate_biodiversity_indices(community_data)

    report_lines.append("Diversity Indices:")
    for index_name, values in diversity_indices.items():
        mean_val = statistics.mean(values) if values else 0
        std_val = statistics.stdev(values) if len(values) > 1 else 0
        report_lines.append(f"  {index_name.capitalize()}: {mean_val:.3f} ± {std_val:.3f}")
    report_lines.append("")

    # Alpha-beta-gamma diversity
    abg = alpha_beta_gamma_diversity(community_data)
    report_lines.append("Alpha-Beta-Gamma Diversity:")
    report_lines.append(f"  Alpha (within-community): {abg['alpha']:.3f}")
    report_lines.append(f"  Beta (between-community): {abg['beta']:.3f}")
    report_lines.append(f"  Gamma (total): {abg['gamma']:.3f}")
    report_lines.append("")

    # Species richness
    richness_values = species_richness(community_data)
    if isinstance(richness_values, list):
        mean_richness = statistics.mean(richness_values)
        report_lines.append(f"Mean Species Richness: {mean_richness:.1f}")

    # Community similarity
    if len(community_data) > 1:
        similarity_matrix = community_similarity_matrix(community_data, "bray_curtis")
        mean_similarity = statistics.mean(
            [
                similarity_matrix[i][j]
                for i in range(len(similarity_matrix))
                for j in range(i + 1, len(similarity_matrix))
            ]
        )
        report_lines.append(f"Mean Community Similarity (Bray-Curtis): {mean_similarity:.3f}")

    report = "\n".join(report_lines)

    if output_path:
        output_path = Path(output_path)
        with open(output_path, "w") as f:
            f.write(report)
        logger.info(f"Ecology report saved to {output_path}")

    return report


def shannon_diversity(abundances: List[float]) -> float:
    """Calculate Shannon diversity index.

    Args:
        abundances: List of species abundances

    Returns:
        Shannon diversity index (H')

    Example:
        >>> shannon_diversity([10, 10, 10])
        1.09861228866811
        >>> shannon_diversity([1, 1, 1, 1])
        1.3862943611198906
    """
    return calculate_single_diversity(abundances, "shannon")


def simpson_diversity(abundances: List[float]) -> float:
    """Calculate Simpson diversity index.

    Args:
        abundances: List of species abundances

    Returns:
        Simpson diversity index (1-D)

    Example:
        >>> simpson_diversity([10, 10, 10])
        0.6666666666666667
        >>> simpson_diversity([1, 1, 1, 1])
        0.75
    """
    return calculate_single_diversity(abundances, "simpson")


def species_richness_simple(abundances: List[float]) -> int:
    """Calculate species richness (count of non-zero abundances).

    Args:
        abundances: List of species abundances

    Returns:
        Number of species with non-zero abundance

    Example:
        >>> species_richness_simple([1, 2, 0, 3])
        3
    """
    if not abundances:
        return 0
    return sum(1 for a in abundances if a > 0)


def pielou_evenness(abundances: List[float]) -> float:
    """Calculate Pielou's evenness index.

    J = H / ln(S) where H is Shannon diversity and S is species richness.

    Args:
        abundances: List of species abundances

    Returns:
        Pielou's evenness index (0 to 1)

    Example:
        >>> pielou_evenness([1, 1, 1, 1])
        1.0
    """
    if not abundances:
        return 0.0

    S = species_richness_simple(abundances)
    if S <= 1:
        return 1.0 if S == 1 else 0.0

    H = shannon_diversity(abundances)
    H_max = math.log(S)

    if H_max == 0:
        return 0.0

    return H / H_max


def chao1_estimator(abundances: List[float]) -> float:
    """Estimate species richness using Chao1 estimator.

    Chao1 = S_obs + f1^2 / (2*f2)  when f2 > 0
    Chao1 = S_obs + f1*(f1-1) / 2  when f2 = 0

    where:
        S_obs = observed species richness
        f1 = number of singletons (species with abundance 1)
        f2 = number of doubletons (species with abundance 2)

    Args:
        abundances: List of species abundances

    Returns:
        Estimated total species richness

    Example:
        >>> chao1_estimator([10, 8, 6, 4, 2, 1, 1, 1])
        12.5
    """
    if not abundances:
        return 0.0

    S_obs = species_richness_simple(abundances)
    if S_obs == 0:
        return 0.0

    f1 = sum(1 for a in abundances if a == 1)  # Singletons
    f2 = sum(1 for a in abundances if a == 2)  # Doubletons

    if f2 > 0:
        return S_obs + (f1**2) / (2 * f2)
    else:
        # Bias-corrected form when f2 = 0
        return S_obs + (f1 * (f1 - 1)) / 2


def community_metrics(abundances: List[float]) -> Dict[str, float]:
    """Calculate comprehensive community ecology metrics.

    Args:
        abundances: List of species abundances

    Returns:
        Dictionary with all community metrics:
            - shannon: Shannon diversity index
            - simpson: Simpson diversity index (1-D)
            - richness: Species richness
            - pielou: Pielou's evenness
            - chao1: Chao1 estimated richness

    Example:
        >>> metrics = community_metrics([10, 8, 6, 4, 2, 1, 1, 1])
        >>> metrics['richness']
        8
    """
    return {
        "shannon": shannon_diversity(abundances),
        "simpson": simpson_diversity(abundances),
        "richness": species_richness_simple(abundances),
        "pielou": pielou_evenness(abundances),
        "chao1": chao1_estimator(abundances),
    }
