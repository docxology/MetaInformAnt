"""Tissue-specificity (tau) and orthology-class stratification methods.

Implements the tau tissue-specificity index (Yanai et al. 2005) following the
protocol of Xu & Colgan 2025 (MBE, 10.1093/molbev/msaf063): TPM-normalized
input, drop of the lowest-10%-mean-expression genes, log2 transform, and
tau = sum(1 - x_i / max(x_i)) / (n - 1) over the n tissues with max > 0.

Also provides pure, table-in/table-out functions to join expression tables
against the ortholog bridge produced by
``projects/hymenoptera_amalgkit/scripts/build_ortholog_bridge.py``
(``transcript_orthogroups.tsv``: rows = orthogroups, columns = species,
cells = comma-separated transcript IDs) and to summarize tau distributions
across orthology classes (single-copy ortholog vs multicopy) descriptively.

Descriptive statistics only in the default outputs. The Wilcoxon
rank-sum comparison of tau across orthology classes is exposed as a separate,
explicitly-gated function for post-evidence-freeze use; it must never be
called on mid-campaign data.
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
import pandas as pd

# =============================================================================
# Tau tissue-specificity
# =============================================================================


def filter_low_expression(
    expression_df: pd.DataFrame,
    lowest_fraction: float = 0.10,
) -> pd.DataFrame:
    """Drop genes in the lowest ``lowest_fraction`` of mean expression.

    Implements the lowest-10%-mean-expression filter from Xu & Colgan 2025.
    Genes whose mean expression equals zero are always excluded. Among
    nonzero genes, the retained set is everything at or above the
    ``lowest_fraction`` quantile of mean expression (ties at the boundary
    kept).

    Args:
        expression_df: Expression matrix (genes x tissues/samples).
        lowest_fraction: Fraction of lowest-mean-expression genes to drop.

    Returns:
        Filtered expression matrix (copy).
    """
    if not 0.0 <= lowest_fraction < 1.0:
        raise ValueError("lowest_fraction must be in [0, 1)")
    if expression_df.empty:
        return expression_df.copy()

    values = expression_df.astype(float)
    mean_expr = values.mean(axis=1)
    nonzero = mean_expr > 0
    if not nonzero.any():
        # Only NaN/zero rows: keep rows containing NaN (NaN discipline) and
        # drop all-zero rows.
        has_nan = values.isna().any(axis=1)
        return values.loc[has_nan].copy()

    kept_mean = mean_expr[nonzero]
    cutoff = kept_mean.quantile(lowest_fraction)
    # Rows containing NaN never enter the ranking; they pass through so tau
    # propagates NaN instead of silently dropping the gene.
    has_nan = values.isna().any(axis=1)
    keep = has_nan | (nonzero & (mean_expr >= cutoff))
    return values.loc[keep].copy()


def compute_tau(
    expression_df: pd.DataFrame,
    log2: bool = True,
    lowest_fraction: float = 0.10,
) -> pd.Series:
    """Compute tau tissue-specificity per gene (Yanai et al. 2005).

    Protocol per Xu & Colgan 2025:
      1. Input is TPM-normalized expression (genes x tissues). If a raw count
         matrix is supplied, normalize first with
         ``metainformant.rna.analysis.expression_core.normalize_counts``.
      2. Drop the lowest-10%-mean-expression genes (and all-zero genes).
      3. Log2-transform (log2(x + 1) to keep zeros finite).
      4. tau = sum_i (1 - x_i / max_i) / (n - 1) over the n tissues where the
         gene's maximum expression is > 0.

    Edge cases:
      - All-zero genes: excluded by the expression filter (tau undefined).
      - Single tissue (n = 1): tau is undefined and returned as NaN.
      - NaN input values: propagate as NaN tau for that gene (no silent
        imputation).

    Args:
        expression_df: TPM expression matrix (genes x tissues).
        log2: Apply the log2(x + 1) transform before computing tau.
        lowest_fraction: Fraction of lowest-mean-expression genes to drop
            before computing tau (default 0.10 per Xu & Colgan 2025).

    Returns:
        Series of tau values indexed by gene. NaN where undefined.
    """
    if expression_df.empty:
        return pd.Series(dtype=float, name="tau")

    filtered = filter_low_expression(expression_df, lowest_fraction=lowest_fraction)
    if filtered.empty:
        return pd.Series(dtype=float, name="tau")

    values = filtered.astype(float)
    if log2:
        values = np.log2(values + 1.0)

    row_max = values.max(axis=1)
    n_tissues = values.shape[1]

    # NaN discipline: NaN rows propagate; do not impute.
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = values.div(row_max.where(row_max > 0), axis=0)
        per_gene = (1.0 - ratios).sum(axis=1, min_count=n_tissues) / (n_tissues - 1)

    tau = per_gene.where(row_max > 0)
    if n_tissues <= 1:
        tau = pd.Series(np.nan, index=values.index, name="tau")
    return tau


def tau_summary(tau: pd.Series) -> Dict[str, object]:
    """Descriptive summary of a tau distribution.

    Args:
        tau: Series of tau values (NaNs allowed and ignored).

    Returns:
        Dict with count, valid_count, nan_count, median, mean, q25, q75.
    """
    valid = tau.dropna()
    return {
        "count": int(tau.size),
        "valid_count": int(valid.size),
        "nan_count": int(tau.isna().sum()),
        "median": float(valid.median()) if not valid.empty else None,
        "mean": float(valid.mean()) if not valid.empty else None,
        "q25": float(valid.quantile(0.25)) if not valid.empty else None,
        "q75": float(valid.quantile(0.75)) if not valid.empty else None,
    }


# =============================================================================
# Orthology-class stratification
# =============================================================================


def classify_orthogroups(
    bridge_table: pd.DataFrame,
    min_species: int = 2,
) -> pd.DataFrame:
    """Classify orthogroups as single-copy orthologs vs multicopy.

    Reads the transcript-orthogroup bridge table produced by
    ``build_ortholog_bridge.py`` (``transcript_orthogroups.tsv``): index =
    orthogroup name, columns = species, cells = comma-separated transcript
    IDs (empty string when unmapped).

    An orthogroup is classified:
      - ``"single_copy"``: mapped in at least ``min_species`` species and
        every mapped species has exactly one transcript.
      - ``"multicopy"``: mapped in at least ``min_species`` species and at
        least one mapped species has more than one transcript.
      - ``"insufficient"``: mapped in fewer than ``min_species`` species.

    Args:
        bridge_table: Orthogroup x species table of comma-separated IDs.
        min_species: Minimum number of mapped species for classification.

    Returns:
        DataFrame indexed by orthogroup with columns ``n_species_mapped``,
        ``max_copies`` (largest per-species transcript count), and
        ``orthology_class``.
    """
    if min_species < 1:
        raise ValueError("min_species must be >= 1")
    if bridge_table.empty:
        return pd.DataFrame(columns=["n_species_mapped", "max_copies", "orthology_class"])

    records = []
    for og_name, row in bridge_table.iterrows():
        n_mapped = 0
        max_copies = 0
        for cell in row:
            ids = _parse_ids(cell)
            if ids:
                n_mapped += 1
                max_copies = max(max_copies, len(ids))
        if n_mapped < min_species:
            cls = "insufficient"
        elif max_copies > 1:
            cls = "multicopy"
        else:
            cls = "single_copy"
        records.append(
            {
                "orthogroup": og_name,
                "n_species_mapped": n_mapped,
                "max_copies": max_copies,
                "orthology_class": cls,
            }
        )
    return pd.DataFrame(records).set_index("orthogroup")


def _parse_ids(cell: object) -> Tuple[str, ...]:
    """Parse a comma-separated bridge-table cell into transcript IDs."""
    if cell is None or (isinstance(cell, float) and np.isnan(cell)):
        return ()
    text = str(cell).strip()
    if not text or text.lower() == "nan":
        return ()
    return tuple(part.strip() for part in text.split(",") if part.strip())


def join_expression_with_orthology(
    expression_df: pd.DataFrame,
    bridge_table: pd.DataFrame,
    species: str,
    min_species: int = 2,
) -> pd.DataFrame:
    """Join an expression table with orthology classes for one species.

    Builds a transcript -> orthogroup map from the bridge table for
    ``species`` and returns a per-transcript frame combining expression
    values with the orthology classification of the transcript's orthogroup.

    Args:
        expression_df: Expression matrix (genes/transcripts x tissues) for
            one species.
        bridge_table: Orthogroup x species bridge table (see
            :func:`classify_orthogroups`).
        species: Column name in ``bridge_table`` matching the expression
            table's species.
        min_species: Minimum mapped species for orthogroup classification.

    Returns:
        DataFrame indexed by transcript with original expression columns
        plus ``orthogroup``, ``orthology_class``, ``n_species_mapped``, and
        ``max_copies``. Transcripts absent from the bridge are labeled
        ``orthology_class`` = ``"unmapped"``.
    """
    if species not in bridge_table.columns:
        raise ValueError(f"Species '{species}' not in bridge table columns")
    classes = classify_orthogroups(bridge_table, min_species=min_species)

    transcript_map: Dict[str, Tuple[str, str]] = {}
    for og_name, cell in bridge_table[species].items():
        for tid in _parse_ids(cell):
            transcript_map[tid] = (str(og_name), str(classes.loc[og_name, "orthology_class"]))

    records = []
    for gene, row in expression_df.iterrows():
        og_name, cls = transcript_map.get(str(gene), ("", "unmapped"))
        record = dict(row)
        record["orthogroup"] = og_name
        record["orthology_class"] = cls
        if og_name:
            info = classes.loc[og_name]
            record["n_species_mapped"] = int(info["n_species_mapped"])
            record["max_copies"] = int(info["max_copies"])
        else:
            record["n_species_mapped"] = 0
            record["max_copies"] = 0
        records.append(record)

    return pd.DataFrame(records, index=expression_df.index)


# =============================================================================
# Duplication-specificity coupling (descriptive only by default)
# =============================================================================


def duplication_specificity_summary(
    tau: pd.Series,
    orthology_class: pd.Series,
) -> pd.DataFrame:
    """Descriptive tau summaries stratified by orthology class.

    Compares tau distributions across orthology classes (e.g.
    single-copy ortholog vs multicopy) with medians and quantiles only.
    No test statistics are produced here; use
    :func:`wilcoxon_duplication_specificity` for the gated inferential
    comparison after the evidence manifest freezes.

    Args:
        tau: Series of tau values indexed by gene.
        orthology_class: Series of class labels, aligned index with ``tau``.

    Returns:
        DataFrame indexed by orthology class with count, valid_count,
        nan_count, median, mean, q25, q75.
    """
    if not tau.index.equals(orthology_class.index):
        raise ValueError("tau and orthology_class must share an index")
    rows = []
    for cls in sorted(orthology_class.dropna().unique()):
        sub = tau[orthology_class == cls]
        rows.append({"orthology_class": cls, **tau_summary(sub)})
    return pd.DataFrame(rows).set_index("orthology_class")


def wilcoxon_duplication_specificity(
    tau: pd.Series,
    orthology_class: pd.Series,
    group_a: str = "single_copy",
    group_b: str = "multicopy",
    evidence_manifest_frozen: bool = False,
) -> Dict[str, object]:
    """GATED: Wilcoxon rank-sum comparison of tau across orthology classes.

    POST-FREEZE USE ONLY. The evidence manifest has not frozen mid-campaign,
    so any inferential claim is out of scope. The
    ``evidence_manifest_frozen`` flag makes the gate explicit at every call
    site: default False always refuses.

    Args:
        tau: Series of tau values indexed by gene.
        orthology_class: Series of class labels aligned with ``tau``.
        group_a: First orthology class (default single_copy).
        group_b: Second orthology class (default multicopy).
        evidence_manifest_frozen: Must be True; the caller affirms the
            evidence manifest has frozen and inferential output is allowed.

    Returns:
        Dict with statistic, p_value, group sizes, medians, and gate label.

    Raises:
        RuntimeError: If ``evidence_manifest_frozen`` is not True or either
            group has fewer than 2 valid observations.
    """
    from scipy.stats import ranksums  # local import: inferential path only

    if not evidence_manifest_frozen:
        raise RuntimeError(
            "wilcoxon_duplication_specificity is gated for post-freeze use: "
            "refusing to run while the evidence manifest is unfrozen"
        )
    if not tau.index.equals(orthology_class.index):
        raise ValueError("tau and orthology_class must share an index")

    a = tau[orthology_class == group_a].dropna()
    b = tau[orthology_class == group_b].dropna()
    if len(a) < 2 or len(b) < 2:
        raise RuntimeError(f"Insufficient observations for Wilcoxon: {group_a}={len(a)}, {group_b}={len(b)}")
    stat, p_value = ranksums(a, b)
    return {
        "statistic": float(stat),
        "p_value": float(p_value),
        "n_" + group_a: int(len(a)),
        "n_" + group_b: int(len(b)),
        "median_" + group_a: float(a.median()),
        "median_" + group_b: float(b.median()),
        "gate": "post-freeze",
    }
