"""Descriptive species mean-profile distances across mapped orthogroups.

This estimand differs from correlating a gene across aligned samples. Mapping
cells select their first recorded transcript, preserving the project's legacy
policy. Unsupported pairs remain NaN, accompanied by shared-orthogroup counts.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def compute_ortholog_profile_divergence(
    transcript_og: pd.DataFrame,
    species_expression: Mapping[str, pd.Series],
    *,
    min_shared_orthologs: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return ``1 - Spearman rho`` and overlap counts without inventing distances."""
    species = list(transcript_og.columns)
    if len(species) < 2 or not transcript_og.columns.is_unique or not transcript_og.index.is_unique:
        raise ValueError("At least two unique mapped species and unique orthogroup IDs are required")
    if type(min_shared_orthologs) is not int or min_shared_orthologs < 2:
        raise ValueError("min_shared_orthologs must be an integer of at least two")
    profiles: dict[str, pd.Series] = {}
    for name in species:
        if name not in species_expression:
            profiles[name] = pd.Series(dtype=float)
            continue
        expression = species_expression[name]
        values = expression.to_numpy(dtype=float)
        if not expression.index.is_unique or not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"Expression profile must have unique IDs and finite nonnegative values: {name}")
        mapped = {}
        for group, row in transcript_og.iterrows():
            cell = row[name]
            if pd.isna(cell) or not str(cell).strip():
                continue
            transcript = str(cell).split(",")[0]
            if transcript in expression.index:
                mapped[group] = expression[transcript]
        profiles[name] = pd.Series(mapped, dtype=float)
    n = len(species)
    distances = np.full((n, n), np.nan)
    np.fill_diagonal(distances, 0)
    counts = np.zeros((n, n), dtype=int)
    for i, left in enumerate(species):
        for j in range(i + 1, n):
            right = species[j]
            shared = profiles[left].index.intersection(profiles[right].index)
            counts[i, j] = counts[j, i] = len(shared)
            if len(shared) < min_shared_orthologs:
                continue
            a, b = profiles[left].loc[shared], profiles[right].loc[shared]
            if a.nunique() < 2 or b.nunique() < 2:
                continue
            rho = float(spearmanr(a, b).statistic)
            if np.isfinite(rho):
                distances[i, j] = distances[j, i] = float(np.clip(1 - rho, 0, 2))
    return (
        pd.DataFrame(distances, index=species, columns=species),
        pd.DataFrame(counts, index=species, columns=species),
    )
