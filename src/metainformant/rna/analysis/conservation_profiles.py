"""Cross-species expression-profile conservation and summary artifacts.

Implements the descriptive, transferable core of the Round-3 methods papers
(Yuan+ 2026 Sci Data; Xu & Colgan 2025 MBE; Leader+ 2024 JMB; Mantica+ 2024
NEE) that the existing modules do not yet provide:

- Spearman-based expression-profile conservation across orthologs
  (:func:`compute_profile_conservation`) — per-gene rank correlation between
  species over matched tissue/condition labels, plus an orthogroup-level
  summary. Descriptive only.
- Cross-species TPM distribution summary artifacts
  (:func:`compute_tpm_distribution_summary`) — per-species/per-tissue
  quantile summaries of TPM-scale expression matrices.
- Per-tissue completeness tables (:func:`compute_per_tissue_completeness`)
  and (:func:`compute_tissue_overlap_summary`) — how many species contribute
  measured expression to each tissue label, the sampling-dependence record
  required before any cross-species tau comparison (claim boundary 2 of
  ``docs/rna/METHODS_LITERATURE_ALIGNMENT.md``).

All outputs are descriptive statistics of the supplied matrices. No p-values,
confidence intervals, or inferential claims are produced here; the evidence
manifest is not frozen.
"""

from __future__ import annotations

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Minimum explicitly matched conditions for a species pair to be scored.
MIN_PROFILE_CONDITIONS = 3

# Quantiles reported by the TPM distribution summary.
TPM_QUANTILES: Tuple[float, ...] = (0.10, 0.25, 0.50, 0.75, 0.90)


# =============================================================================
# Expression-profile conservation across orthologs (descriptive)
# =============================================================================


def _validate_species_profiles(
    species_profiles: Mapping[str, pd.DataFrame],
    min_species: int = 2,
) -> None:
    """Reject malformed species profile mappings before analysis."""

    if len(species_profiles) < min_species:
        raise ValueError(f"analysis requires at least {min_species} species")
    for species, profile in species_profiles.items():
        if not isinstance(profile, pd.DataFrame):
            raise TypeError(f"profile for species '{species}' must be a pandas DataFrame")
        if profile.columns.duplicated().any():
            raise ValueError(f"profile for species '{species}' has duplicate condition labels")
        if profile.index.duplicated().any():
            raise ValueError(f"profile for species '{species}' has duplicate gene labels")


def compute_profile_conservation(
    species_profiles: Mapping[str, pd.DataFrame],
    profile_alignments: Optional[Mapping[str, str]] = None,
    method: str = "spearman",
) -> pd.DataFrame:
    """Compute pairwise expression-profile conservation across species.

    For every species pair, aligns the condition/tissue columns with
    ``profile_alignments`` (species condition label -> canonical label) and
    computes, for each gene present in both species, the correlation of its
    expression profile across the matched conditions. This is the
    Mantica+ 2024 / Xu & Colgan 2025 style per-gene profile comparison,
    restricted to explicitly matched tissue labels: no positional alignment
    is ever inferred.

    Args:
        species_profiles: Mapping species -> genes x conditions DataFrame
            (TPM-scale or any comparable per-gene scale).
        profile_alignments: Optional mapping from each species' condition
            labels to a shared canonical label. Labels that already match may
            be omitted. Required when species use different labels.
        method: "spearman" (default, rank-based — the papers' convention) or
            "pearson".

    Returns:
        Long-format DataFrame with columns: species_a, species_b, gene_id,
        n_conditions, correlation. Pairs with fewer than
        ``MIN_PROFILE_CONDITIONS`` matched conditions, genes absent from one
        species, constant profiles, and non-computable correlations are
        dropped deterministically. Sorted by species pair then gene.

    Raises:
        ValueError: If fewer than two species are supplied, a profile is
            malformed, or ``method`` is unknown.
    """

    if method not in ("spearman", "pearson"):
        raise ValueError(f"Unknown method: {method}. Valid: spearman, pearson")
    _validate_species_profiles(species_profiles)

    species_names = sorted(species_profiles)
    align = dict(profile_alignments or {})

    canonical_sets: Dict[str, set] = {}
    inverse_maps: Dict[str, Dict[str, str]] = {}
    for species in species_names:
        inverse: Dict[str, str] = {}
        for label in species_profiles[species].columns:
            canonical = str(align.get(label, label))
            if canonical in inverse:
                raise ValueError(
                    f"species '{species}' maps condition labels "
                    f"'{inverse[canonical]}' and '{label}' to the same "
                    f"canonical label '{canonical}'"
                )
            inverse[canonical] = str(label)
        inverse_maps[species] = inverse
        canonical_sets[species] = set(inverse)

    records: List[Dict[str, object]] = []
    for i, sp_a in enumerate(species_names):
        for sp_b in species_names[i + 1 :]:
            shared = sorted(canonical_sets[sp_a] & canonical_sets[sp_b])
            if len(shared) < MIN_PROFILE_CONDITIONS:
                logger.warning(
                    "Pair (%s, %s): only %d matched conditions; skipping",
                    sp_a,
                    sp_b,
                    len(shared),
                )
                continue

            a_labels = [inverse_maps[sp_a][c] for c in shared]
            b_labels = [inverse_maps[sp_b][c] for c in shared]

            sub_a = species_profiles[sp_a][a_labels]
            sub_b = species_profiles[sp_b][b_labels]
            shared_genes = sorted(set(sub_a.index) & set(sub_b.index))
            if not shared_genes:
                logger.warning("Pair (%s, %s): no shared genes; skipping", sp_a, sp_b)
                continue

            vec_a = sub_a.loc[shared_genes].to_numpy(dtype=float)
            vec_b = sub_b.loc[shared_genes].to_numpy(dtype=float)

            corr_func = stats.spearmanr if method == "spearman" else stats.pearsonr
            for row_idx, gene in enumerate(shared_genes):
                x = vec_a[row_idx]
                y = vec_b[row_idx]
                if np.std(x) == 0.0 or np.std(y) == 0.0:
                    continue
                result = corr_func(x, y)
                correlation = float(result.statistic if hasattr(result, "statistic") else result[0])
                if np.isnan(correlation):
                    continue
                records.append(
                    {
                        "species_a": sp_a,
                        "species_b": sp_b,
                        "gene_id": gene,
                        "n_conditions": len(shared),
                        "correlation": correlation,
                    }
                )

    result_df = pd.DataFrame(
        records,
        columns=["species_a", "species_b", "gene_id", "n_conditions", "correlation"],
    )
    if result_df.empty:
        return result_df
    return result_df.sort_values(["species_a", "species_b", "gene_id"], kind="mergesort").reset_index(drop=True)


def summarize_profile_conservation(profile_conservation: pd.DataFrame) -> pd.DataFrame:
    """Aggregate long-format profile conservation to one row per gene.

    Descriptive summary: count, mean, median, min, and max of the pairwise
    profile correlations for each gene across the species pairs in which it
    could be scored. No significance is attached to any value.

    Args:
        profile_conservation: Output of :func:`compute_profile_conservation`.

    Returns:
        DataFrame indexed by ``gene_id`` with columns n_pairs,
        mean_correlation, median_correlation, min_correlation,
        max_correlation, sorted by mean descending. Empty input yields the
        same columns with zero rows.
    """

    columns = [
        "n_pairs",
        "mean_correlation",
        "median_correlation",
        "min_correlation",
        "max_correlation",
    ]
    if profile_conservation.empty:
        return pd.DataFrame(columns=columns, index=pd.Index([], name="gene_id"))

    grouped = profile_conservation.groupby("gene_id")["correlation"]
    summary = pd.DataFrame(
        {
            "n_pairs": grouped.count(),
            "mean_correlation": grouped.mean(),
            "median_correlation": grouped.median(),
            "min_correlation": grouped.min(),
            "max_correlation": grouped.max(),
        }
    )
    summary.index.name = "gene_id"
    return summary.sort_values("mean_correlation", ascending=False)


# =============================================================================
# Cross-species TPM distribution summaries
# =============================================================================


def compute_tpm_distribution_summary(
    species_profiles: Mapping[str, pd.DataFrame],
    quantiles: Sequence[float] = TPM_QUANTILES,
) -> pd.DataFrame:
    """Summarize the TPM distribution per species and per tissue.

    Descriptive quantile summary of every species x condition cell in the
    supplied profiles. This is the distribution-level artifact the atlas
    grammar (Leader+ 2024) and any post-freeze distributional comparison
    require as their documented baseline.

    Args:
        species_profiles: Mapping species -> genes x conditions DataFrame of
            TPM-scale values.
        quantiles: Quantiles to report (default 10/25/50/75/90%).

    Returns:
        Long-format DataFrame with columns: species, condition, n_genes,
        n_expressed, mean_tpm, median_tpm, plus one ``q{int(q*100)}`` column
        per requested quantile. Sorted deterministically by species and
        condition.

    Raises:
        ValueError: If a quantile falls outside (0, 1) or a profile is not a
            DataFrame.
    """

    for q in quantiles:
        if not 0.0 < q < 1.0:
            raise ValueError(f"quantiles must be in (0, 1); got {q}")
    _validate_species_profiles(species_profiles, min_species=1)

    records: List[Dict[str, object]] = []
    for species in sorted(species_profiles):
        profile = species_profiles[species].astype(float)
        for condition in sorted(profile.columns, key=str):
            values = profile[condition].to_numpy(dtype=float)
            finite = values[np.isfinite(values)]
            n_genes = int(finite.size)
            n_expressed = int((finite > 0).sum())
            record: Dict[str, object] = {
                "species": species,
                "condition": str(condition),
                "n_genes": n_genes,
                "n_expressed": n_expressed,
                "mean_tpm": round(float(finite.mean()), 4) if n_genes else np.nan,
                "median_tpm": round(float(np.median(finite)), 4) if n_genes else np.nan,
            }
            for q in quantiles:
                record[f"q{int(round(q * 100))}"] = round(float(np.quantile(finite, q)), 4) if n_genes else np.nan
            records.append(record)
    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records).sort_values(["species", "condition"], kind="mergesort").reset_index(drop=True)


# =============================================================================
# Per-tissue completeness
# =============================================================================


def compute_per_tissue_completeness(
    species_profiles: Mapping[str, pd.DataFrame],
    profile_alignments: Optional[Mapping[str, str]] = None,
    min_expression: float = 0.0,
) -> pd.DataFrame:
    """Build the per-tissue completeness table across species.

    For each canonical tissue label, reports which species contribute
    measured expression (any gene above ``min_expression``) and the count of
    measured genes per species. This is the sampling-dependence record
    demanded by the claim-boundary rules: cross-species tau comparisons are
    admissible only for tissue sets whose completeness this table documents.

    Args:
        species_profiles: Mapping species -> genes x conditions DataFrame.
        profile_alignments: Optional species-condition -> canonical-tissue
            mapping (same convention as :func:`compute_profile_conservation`).
        min_expression: Minimum value in a tissue column for the tissue to
            count as measured for that species (default: any nonzero value).

    Returns:
        DataFrame indexed by canonical tissue with one boolean column per
        species (``measured_<species>``), one integer column per species
        (``n_genes_<species>``), and ``n_species_measured``. Sorted by tissue.
    """

    _validate_species_profiles(species_profiles)
    align = dict(profile_alignments or {})
    species_names = sorted(species_profiles)

    measured: Dict[str, Dict[str, bool]] = {}
    gene_counts: Dict[str, Dict[str, int]] = {}
    for species in species_names:
        profile = species_profiles[species].astype(float)
        for condition in profile.columns:
            tissue = str(align.get(condition, condition))
            column = profile[condition]
            n_measured = int((column > min_expression).sum())
            measured.setdefault(tissue, {})[species] = n_measured > 0
            gene_counts.setdefault(tissue, {})[species] = max(
                gene_counts.setdefault(tissue, {}).get(species, 0), n_measured
            )

    records: List[Dict[str, object]] = []
    for tissue in sorted(measured):
        record: Dict[str, object] = {"tissue": tissue}
        for species in species_names:
            record[f"measured_{species}"] = bool(measured[tissue].get(species, False))
            record[f"n_genes_{species}"] = int(gene_counts[tissue].get(species, 0))
        record["n_species_measured"] = sum(bool(measured[tissue].get(species, False)) for species in species_names)
        records.append(record)
    if not records:
        return pd.DataFrame().set_index(pd.Index([], name="tissue"))
    return pd.DataFrame(records).set_index("tissue")


def compute_tissue_overlap_summary(
    species_profiles: Mapping[str, pd.DataFrame],
    profile_alignments: Optional[Mapping[str, str]] = None,
) -> pd.DataFrame:
    """Summarize matched-tissue coverage per species pair.

    Complements :func:`compute_per_tissue_completeness` with the pairwise
    view the profile-conservation module needs: for each species pair, how
    many canonical tissue labels are matched on both sides. Pairs below
    ``MIN_PROFILE_CONDITIONS`` matched tissues are exactly the ones
    :func:`compute_profile_conservation` will skip, so this table lets
    callers predict that deterministically.

    Args:
        species_profiles: Mapping species -> genes x conditions DataFrame.
        profile_alignments: Optional canonical-tissue mapping.

    Returns:
        DataFrame with columns: species_a, species_b, n_matched_tissues,
        matched_tissues (comma-separated, deterministic order). One row per
        unordered pair, sorted by species names.
    """

    _validate_species_profiles(species_profiles)
    align = dict(profile_alignments or {})
    species_names = sorted(species_profiles)

    tissue_sets: Dict[str, set] = {}
    for species in species_names:
        tissue_sets[species] = {str(align.get(condition, condition)) for condition in species_profiles[species].columns}

    records: List[Dict[str, object]] = []
    for i, sp_a in enumerate(species_names):
        for sp_b in species_names[i + 1 :]:
            matched = sorted(tissue_sets[sp_a] & tissue_sets[sp_b])
            records.append(
                {
                    "species_a": sp_a,
                    "species_b": sp_b,
                    "n_matched_tissues": len(matched),
                    "matched_tissues": ",".join(matched),
                }
            )
    return pd.DataFrame(records, columns=["species_a", "species_b", "n_matched_tissues", "matched_tissues"])


__all__ = [
    "MIN_PROFILE_CONDITIONS",
    "TPM_QUANTILES",
    "compute_per_tissue_completeness",
    "compute_profile_conservation",
    "compute_tissue_overlap_summary",
    "compute_tpm_distribution_summary",
    "summarize_profile_conservation",
]
