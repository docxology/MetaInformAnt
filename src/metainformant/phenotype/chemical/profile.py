"""Chemical profile analysis for GC-MS, CHC, and metabolomic data."""

from __future__ import annotations

import math
from typing import Dict, List, Any, Optional, Tuple
from collections import defaultdict

from .compound import Compound
from metainformant.core.utils.errors import ValidationError


class ChemicalProfile:
    """Represents a chemical profile (e.g., from GC-MS).

    Attributes:
        sample_id: Unique identifier for the sample.
        compounds: Mapping of Compound to abundance/concentration.
        metadata: Experimental conditions.
    """

    def __init__(
        self, sample_id: str, compounds: Dict[Compound, float], metadata: Dict[str, Any] | None = None
    ) -> None:
        self.sample_id = sample_id
        self.compounds = compounds
        self.metadata = metadata or {}

    @property
    def n_compounds(self) -> int:
        return len(self.compounds)

    @property
    def total_abundance(self) -> float:
        return sum(self.compounds.values())

    def normalize(self, method: str = "total_ion_current") -> ChemicalProfile:
        """Return a new normalized profile.

        Args:
            method: 'total_ion_current' (sum to 1) or 'max_peak' (max to 1).
        """
        if not self.compounds:
            return ChemicalProfile(self.sample_id, {}, self.metadata)

        if method == "total_ion_current":
            total = sum(self.compounds.values())
            if total == 0:
                return ChemicalProfile(self.sample_id, self.compounds.copy(), self.metadata)
            new_compounds = {k: v / total for k, v in self.compounds.items()}
        elif method == "max_peak":
            max_val = max(self.compounds.values())
            if max_val == 0:
                return ChemicalProfile(self.sample_id, self.compounds.copy(), self.metadata)
            new_compounds = {k: v / max_val for k, v in self.compounds.items()}
        else:
            raise ValidationError(f"Unknown normalization method: {method}")

        return ChemicalProfile(self.sample_id, new_compounds, self.metadata)

    def bray_curtis_distance(self, other: ChemicalProfile) -> float:
        """Bray-Curtis dissimilarity between two profiles. Range [0, 1]."""
        all_compounds = set(self.compounds.keys()) | set(other.compounds.keys())
        diff_sum = 0.0
        total_sum = 0.0

        for cmp in all_compounds:
            val1 = self.compounds.get(cmp, 0.0)
            val2 = other.compounds.get(cmp, 0.0)
            diff_sum += abs(val1 - val2)
            total_sum += val1 + val2

        if total_sum == 0:
            return 0.0
        return diff_sum / total_sum

    def euclidean_distance(self, other: ChemicalProfile) -> float:
        """Euclidean distance between two profiles."""
        all_compounds = set(self.compounds.keys()) | set(other.compounds.keys())
        sq_sum = 0.0
        for cmp in all_compounds:
            diff = self.compounds.get(cmp, 0.0) - other.compounds.get(cmp, 0.0)
            sq_sum += diff * diff
        return math.sqrt(sq_sum)

    def cosine_similarity(self, other: ChemicalProfile) -> float:
        """Cosine similarity between two profiles. Range [-1, 1]."""
        all_compounds = set(self.compounds.keys()) | set(other.compounds.keys())
        dot = 0.0
        norm_a = 0.0
        norm_b = 0.0

        for cmp in all_compounds:
            a = self.compounds.get(cmp, 0.0)
            b = other.compounds.get(cmp, 0.0)
            dot += a * b
            norm_a += a * a
            norm_b += b * b

        denom = math.sqrt(norm_a) * math.sqrt(norm_b)
        if denom == 0:
            return 0.0
        return dot / denom

    def shannon_diversity(self) -> float:
        """Shannon diversity index of compound abundances."""
        total = sum(self.compounds.values())
        if total == 0:
            return 0.0

        h = 0.0
        for v in self.compounds.values():
            p = v / total
            if p > 0:
                h -= p * math.log(p)
        return h

    def richness(self, threshold: float = 0.0) -> int:
        """Number of compounds above an abundance threshold."""
        return sum(1 for v in self.compounds.values() if v > threshold)

    def top_compounds(self, n: int = 10) -> List[Tuple[Compound, float]]:
        """Return top-n compounds by abundance."""
        sorted_compounds = sorted(self.compounds.items(), key=lambda x: x[1], reverse=True)
        return sorted_compounds[:n]

    def filter_by_abundance(self, min_abundance: float = 0.0, max_abundance: float = float("inf")) -> ChemicalProfile:
        """Return a new profile with compounds within abundance range."""
        filtered = {c: v for c, v in self.compounds.items() if min_abundance <= v <= max_abundance}
        return ChemicalProfile(self.sample_id, filtered, self.metadata)

    def relative_abundance(self, compound: Compound) -> float:
        """Relative abundance of a specific compound (fraction of total)."""
        total = sum(self.compounds.values())
        if total == 0:
            return 0.0
        return self.compounds.get(compound, 0.0) / total

    def compound_ratio(self, numerator: Compound, denominator: Compound) -> Optional[float]:
        """Ratio between two compounds. Returns None if denominator is absent or zero."""
        num_val = self.compounds.get(numerator, 0.0)
        den_val = self.compounds.get(denominator, 0.0)
        if den_val == 0:
            return None
        return num_val / den_val


def distance_matrix(profiles: List[ChemicalProfile], metric: str = "bray_curtis") -> List[List[float]]:
    """Compute pairwise distance matrix for a list of chemical profiles.

    Args:
        profiles: List of ChemicalProfile objects.
        metric: 'bray_curtis', 'euclidean', or 'cosine'.

    Returns:
        Symmetric distance matrix as list of lists.
    """
    n = len(profiles)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            if metric == "bray_curtis":
                d = profiles[i].bray_curtis_distance(profiles[j])
            elif metric == "euclidean":
                d = profiles[i].euclidean_distance(profiles[j])
            elif metric == "cosine":
                d = 1.0 - profiles[i].cosine_similarity(profiles[j])
            else:
                raise ValidationError(f"Unknown metric: {metric}")
            matrix[i][j] = d
            matrix[j][i] = d

    return matrix


def identify_marker_compounds(
    group_a: List[ChemicalProfile],
    group_b: List[ChemicalProfile],
    min_fold_change: float = 2.0,
) -> List[Dict[str, Any]]:
    """Identify compounds that differ significantly between two groups.

    Uses mean abundance fold-change to flag differential compounds.

    Args:
        group_a: Profiles from group A.
        group_b: Profiles from group B.
        min_fold_change: Minimum fold-change to report.

    Returns:
        List of dicts with compound, mean_a, mean_b, fold_change.
    """
    all_compounds: set[Compound] = set()
    for p in group_a + group_b:
        all_compounds.update(p.compounds.keys())

    markers = []
    for compound in all_compounds:
        vals_a = [p.compounds.get(compound, 0.0) for p in group_a]
        vals_b = [p.compounds.get(compound, 0.0) for p in group_b]

        mean_a = sum(vals_a) / len(vals_a) if vals_a else 0.0
        mean_b = sum(vals_b) / len(vals_b) if vals_b else 0.0

        denom = max(mean_a, mean_b)
        numer = min(mean_a, mean_b)

        if numer > 0:
            fc = denom / numer
        elif denom > 0:
            fc = float("inf")
        else:
            continue

        if fc >= min_fold_change:
            markers.append(
                {
                    "compound": compound,
                    "mean_a": mean_a,
                    "mean_b": mean_b,
                    "fold_change": fc,
                    "enriched_in": "A" if mean_a > mean_b else "B",
                }
            )

    markers.sort(key=lambda m: m["fold_change"], reverse=True)
    return markers
