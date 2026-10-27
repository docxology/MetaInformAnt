"""Interactive visualization functions for GWAS using Plotly.

This module provides interactive HTML-based plots for GWAS results, including
Manhattan plots, PCA 3D scatter plots, and volcano plots. All functions
gracefully fall back to producing styled HTML tables when Plotly is unavailable.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    import plotly.express as px
    import plotly.graph_objects as go

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    go = None
    px = None


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _ensure_output_dir(output_file: Optional[Union[str, Path]]) -> Optional[Path]:
    """Resolve output file path and create parent directories."""
    if output_file is None:
        return None
    out = Path(output_file)
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def _fallback_html_table(
    title: str,
    headers: List[str],
    rows: List[List[str]],
    output_file: Path,
    max_rows: int = 50,
) -> None:
    """Write a self-contained HTML page with a styled table of top results.

    This is used as a fallback when Plotly is not installed.
    """
    truncated_rows = rows[:max_rows]

    header_cells = "".join(f"<th>{h}</th>" for h in headers)
    body_rows = ""
    for row in truncated_rows:
        cells = "".join(f"<td>{c}</td>" for c in row)
        body_rows += f"<tr>{cells}</tr>\n"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{title}</title>
<style>
body {{ font-family: Arial, Helvetica, sans-serif; margin: 2rem; background: #fafafa; }}
h1 {{ color: #333; }}
table {{ border-collapse: collapse; width: 100%; margin-top: 1rem; }}
th, td {{ padding: 8px 12px; border: 1px solid #ddd; text-align: left; }}
th {{ background: #4a90d9; color: white; }}
tr:nth-child(even) {{ background: #f2f2f2; }}
tr:hover {{ background: #e0ecf8; }}
.note {{ color: #888; font-size: 0.9rem; margin-top: 0.5rem; }}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="note">Showing top {len(truncated_rows)} of {len(rows)} entries.
Install <code>plotly</code> for interactive plots.</p>
<table>
<thead><tr>{header_cells}</tr></thead>
<tbody>
{body_rows}
</tbody>
</table>
</body>
</html>"""

    output_file.write_text(html, encoding="utf-8")


# ---------------------------------------------------------------------------
# 1. Interactive Manhattan Plot
# ---------------------------------------------------------------------------

#: Result keys that may hold the chromosome, in priority order (``chrom`` is what
#: ``run_*_gwas`` results and VCF variant records use).
_CHROMOSOME_KEYS = ("chromosome", "chrom", "chr", "CHROM", "CHR")
_POSITION_KEYS = ("position", "pos", "POS")
_P_VALUE_KEYS = ("p_value", "pvalue", "P")


def _first_present(record: Dict[str, Any], keys: Tuple[str, ...]) -> Any:
    """Return the first value among ``keys`` that is present and not None."""
    for key in keys:
        value = record.get(key)
        if value is not None:
            return value
    return None


def _chromosome_label(value: Any) -> Optional[str]:
    """Normalise a chromosome value to a label; None if absent or blank.

    Integral floats (``1.0``, as read from a numeric column with gaps) become
    ``"1"`` so they sort and group with ``"1"``.
    """
    if value is None:
        return None
    if isinstance(value, float):
        if not math.isfinite(value):
            return None
        if value.is_integer():
            value = int(value)
    label = str(value).strip()
    return label or None


def _integer_position(value: Any) -> Optional[int]:
    """Convert a base-pair position to ``int``; None if absent or not a finite number."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if math.isfinite(number) else None


def _positive_p_value(value: Any) -> Optional[float]:
    """Convert a p-value to ``float``; None unless it is finite and positive."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and number > 0 else None


def _chrom_sort_key(c: str) -> Tuple[int, str]:
    """Sort numeric chromosomes numerically, then the rest (X, Y, MT, ...) by name."""
    c_clean = c.replace("chr", "")
    try:
        return (0, str(int(c_clean)).zfill(5))
    except ValueError:
        return (1, c_clean)


def _prepare_manhattan_records(
    assoc_results: List[Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], List[str], int]:
    """Normalise association results into plottable Manhattan records.

    A result is usable only if it has a chromosome, a position and a positive
    finite p-value; nothing is defaulted, so a variant is never drawn on a
    chromosome it was not assigned to. Each record gets its cumulative x position
    (``cum_pos``, chromosomes laid end to end in sorted order) and ``neg_log_p``.

    Args:
        assoc_results: Association result dicts (see :func:`interactive_manhattan`).

    Returns:
        Tuple of (records, chromosome labels in plotting order, number of results
        skipped as unusable).
    """
    records: List[Dict[str, Any]] = []
    n_skipped = 0
    for r in assoc_results:
        chrom = _chromosome_label(_first_present(r, _CHROMOSOME_KEYS))
        pos = _integer_position(_first_present(r, _POSITION_KEYS))
        pval = _positive_p_value(_first_present(r, _P_VALUE_KEYS))
        if chrom is None or pos is None or pval is None:
            n_skipped += 1
            continue
        vid = _first_present(r, ("variant_id", "rsid"))
        if vid is None:
            vid = f"{chrom}:{pos}"
        records.append(
            {
                "chromosome": chrom,
                "position": pos,
                "p_value": pval,
                "variant_id": str(vid),
                "beta": _first_present(r, ("beta", "effect_size")),
            }
        )

    chrom_set = sorted({r["chromosome"] for r in records}, key=_chrom_sort_key)

    # Per-chromosome max position gives the cumulative offset of the next chromosome
    chrom_max: Dict[str, int] = {}
    for r in records:
        c = r["chromosome"]
        chrom_max[c] = max(chrom_max.get(c, 0), r["position"])

    cum_offset: Dict[str, int] = {}
    running = 0
    for c in chrom_set:
        cum_offset[c] = running
        running += chrom_max.get(c, 0) + 1

    for r in records:
        r["cum_pos"] = r["position"] + cum_offset[r["chromosome"]]
        r["neg_log_p"] = -math.log10(max(r["p_value"], 1e-300))

    return records, chrom_set, n_skipped


def interactive_manhattan(
    assoc_results: List[Dict[str, Any]],
    output_file: Optional[Union[str, Path]] = None,
    gene_annotations: Optional[List[Dict[str, Any]]] = None,
    significance_threshold: float = 5e-8,
    title: str = "Interactive Manhattan Plot",
) -> Dict[str, Any]:
    """Create an interactive Manhattan plot of GWAS association results.

    Generates a Plotly scatter of -log10(p) by cumulative genomic position with
    alternate chromosome colouring, hover information, and a significance
    threshold line.  Falls back to a styled HTML table of top hits when Plotly
    is not installed.

    Args:
        assoc_results: List of dicts, each with at least a chromosome
            (``"chromosome"``, ``"chrom"``, ``"chr"``, ``"CHROM"`` or ``"CHR"``), a
            position (``"position"``, ``"pos"`` or ``"POS"``), and a p-value
            (``"p_value"``, ``"pvalue"`` or ``"P"``).  Optional keys:
            ``"variant_id"``, ``"beta"``.  Results missing any of the three
            required values are skipped (and counted in ``"n_skipped"``), never
            placed on a default chromosome.
        output_file: Path for the output HTML file.
        gene_annotations: Optional list of dicts with ``"variant_id"`` and
            ``"gene"`` keys to annotate specific points.
        significance_threshold: Horizontal line for genome-wide significance.
        title: Plot title.

    Returns:
        Dict with ``"status"``, ``"output_path"``, ``"n_variants"``,
        ``"n_skipped"``, and ``"interactive"`` keys.
    """
    if not assoc_results:
        logger.warning("No association results provided for Manhattan plot")
        return {
            "status": "skipped",
            "output_path": None,
            "n_variants": 0,
            "n_skipped": 0,
            "interactive": False,
        }

    out = _ensure_output_dir(output_file)

    # --- Normalise keys, sort chromosomes and compute cumulative positions
    records, chrom_set, n_skipped = _prepare_manhattan_records(assoc_results)
    if n_skipped:
        logger.warning(
            f"Skipped {n_skipped} of {len(assoc_results)} results without a chromosome, "
            "position and positive p-value in the Manhattan plot"
        )

    if not records:
        logger.warning("No valid variants after filtering for Manhattan plot")
        return {
            "status": "skipped",
            "output_path": None,
            "n_variants": 0,
            "n_skipped": n_skipped,
            "interactive": False,
        }

    n_variants = len(records)
    chrom_to_idx = {c: i for i, c in enumerate(chrom_set)}

    # --- Build annotations lookup if provided --------------------------
    annotation_map: Dict[str, str] = {}
    if gene_annotations:
        for ga in gene_annotations:
            vid = ga.get("variant_id", "")
            gene = ga.get("gene", "")
            if vid and gene:
                annotation_map[str(vid)] = gene

    # --- Plotly path ---------------------------------------------------
    if HAS_PLOTLY and go is not None:
        fig = go.Figure()

        # Discrete colour palette for chromosomes
        palette = px.colors.qualitative.D3 if px is not None else ["#1f77b4", "#ff7f0e"]

        for c in chrom_set:
            chrom_records = [r for r in records if r["chromosome"] == c]
            colour = palette[chrom_to_idx[c] % len(palette)]

            xs = [r["cum_pos"] for r in chrom_records]
            ys = [r["neg_log_p"] for r in chrom_records]
            hover_texts = []
            for r in chrom_records:
                parts = [
                    f"Variant: {r['variant_id']}",
                    f"Chr: {r['chromosome']}",
                    f"Pos: {r['position']:,}",
                    f"P-value: {r['p_value']:.2e}",
                ]
                if r["beta"] is not None:
                    parts.append(f"Beta: {r['beta']:.4f}")
                hover_texts.append("<br>".join(parts))

            fig.add_trace(
                go.Scattergl(
                    x=xs,
                    y=ys,
                    mode="markers",
                    marker=dict(color=colour, size=4, opacity=0.7),
                    text=hover_texts,
                    hoverinfo="text",
                    name=f"Chr {c}",
                    showlegend=False,
                )
            )

        # Significance threshold line
        sig_line_y = -math.log10(significance_threshold)
        fig.add_hline(
            y=sig_line_y,
            line_dash="dash",
            line_color="red",
            annotation_text=f"p = {significance_threshold:.0e}",
        )

        # Gene annotations
        if annotation_map:
            for r in records:
                gene = annotation_map.get(r["variant_id"])
                if gene:
                    fig.add_annotation(
                        x=r["cum_pos"],
                        y=r["neg_log_p"],
                        text=gene,
                        showarrow=True,
                        arrowhead=2,
                        ax=0,
                        ay=-30,
                        font=dict(size=10, color="darkred"),
                    )

        # Axis labels with chromosome tick marks
        chrom_centers = []
        for c in chrom_set:
            chrom_records_pos = [r["cum_pos"] for r in records if r["chromosome"] == c]
            if chrom_records_pos:
                chrom_centers.append((c, (min(chrom_records_pos) + max(chrom_records_pos)) / 2))

        fig.update_layout(
            title=title,
            xaxis=dict(
                title="Chromosome",
                tickvals=[cc[1] for cc in chrom_centers],
                ticktext=[cc[0] for cc in chrom_centers],
            ),
            yaxis=dict(title="-log<sub>10</sub>(p-value)"),
            hovermode="closest",
            template="plotly_white",
        )

        if out is not None:
            fig.write_html(str(out), include_plotlyjs=True)
            logger.info(f"Saved interactive Manhattan plot to {out}")

        return {
            "status": "success",
            "output_path": str(out) if out else None,
            "n_variants": n_variants,
            "n_skipped": n_skipped,
            "interactive": True,
        }

    # --- Fallback: styled HTML table -----------------------------------
    logger.info("Plotly not available; generating fallback HTML table for Manhattan plot")
    sorted_records = sorted(records, key=lambda r: r["p_value"])
    headers = ["Variant", "Chromosome", "Position", "P-value", "Beta", "-log10(p)"]
    rows = []
    for r in sorted_records:
        rows.append(
            [
                r["variant_id"],
                r["chromosome"],
                str(r["position"]),
                f"{r['p_value']:.2e}",
                f"{r['beta']:.4f}" if r["beta"] is not None else "N/A",
                f"{r['neg_log_p']:.2f}",
            ]
        )
    if out is not None:
        _fallback_html_table(title, headers, rows, out)
        logger.info(f"Saved fallback Manhattan table to {out}")

    return {
        "status": "success",
        "output_path": str(out) if out else None,
        "n_variants": n_variants,
        "n_skipped": n_skipped,
        "interactive": False,
    }


# ---------------------------------------------------------------------------
# 2. Interactive 3D PCA
# ---------------------------------------------------------------------------


def interactive_pca(
    pca_data: Dict[str, Any],
    metadata: Optional[Dict[str, Dict]] = None,
    output_file: Optional[Union[str, Path]] = None,
    color_by: str = "population",
    title: str = "Interactive 3D PCA",
) -> Dict[str, Any]:
    """Create an interactive 3D PCA scatter plot.

    Uses Plotly to render PC1, PC2, PC3 in a rotatable 3D scatter. Points are
    coloured by a metadata field (default ``"population"``).  Falls back to a
    styled HTML summary table when Plotly is not installed.

    Args:
        pca_data: Dictionary with:
            - ``"pcs"``: list of lists or 2-D array (samples x components).
            - ``"explained_variance_ratio"``: list of floats per component.
            - ``"sample_ids"`` (optional): list of sample identifiers.
        metadata: Optional mapping ``{sample_id: {field: value, ...}}``.
        output_file: Path for the output HTML file.
        color_by: Metadata field name used for colouring points.
        title: Plot title.

    Returns:
        Dict with ``"status"``, ``"output_path"``, ``"n_samples"``, and
        ``"interactive"`` keys.
    """
    pcs_raw = pca_data.get("pcs", [])
    explained = pca_data.get("explained_variance_ratio", [])
    sample_ids = pca_data.get("sample_ids", [])

    if not pcs_raw or (hasattr(pcs_raw, "__len__") and len(pcs_raw) == 0):
        logger.warning("No PCA data provided")
        return {
            "status": "skipped",
            "output_path": None,
            "n_samples": 0,
            "interactive": False,
        }

    # Convert to list-of-lists for uniform handling
    try:
        import numpy as np

        pcs_arr = np.asarray(pcs_raw, dtype=np.float64)
    except Exception:
        pcs_arr = None

    if pcs_arr is None or pcs_arr.ndim != 2:
        logger.error("pca_data['pcs'] must be a 2D structure")
        return {
            "status": "failed",
            "output_path": None,
            "n_samples": 0,
            "interactive": False,
        }

    n_samples, n_pcs = pcs_arr.shape
    if n_pcs < 3:
        logger.error(f"Need at least 3 PCs for 3D PCA; got {n_pcs}")
        return {
            "status": "failed",
            "output_path": None,
            "n_samples": n_samples,
            "interactive": False,
        }

    out = _ensure_output_dir(output_file)

    # Resolve per-sample labels
    if not sample_ids:
        sample_ids = [f"Sample_{i}" for i in range(n_samples)]

    pop_labels: List[str] = []
    if metadata:
        for sid in sample_ids:
            entry = metadata.get(sid, {})
            pop_labels.append(str(entry.get(color_by, "unknown")))
    else:
        pop_labels = ["all"] * n_samples

    # Axis labels with explained variance
    def _axis_label(pc_idx: int) -> str:
        base = f"PC{pc_idx + 1}"
        if pc_idx < len(explained):
            base += f" ({explained[pc_idx] * 100:.1f}%)"
        return base

    # --- Plotly path ---------------------------------------------------
    if HAS_PLOTLY and go is not None:
        fig = go.Figure()

        unique_pops = sorted(set(pop_labels))
        palette = px.colors.qualitative.D3 if px is not None else ["#1f77b4", "#ff7f0e"]

        for i, pop in enumerate(unique_pops):
            mask = [p == pop for p in pop_labels]
            indices = [j for j, m in enumerate(mask) if m]
            hover_texts = []
            for j in indices:
                parts = [
                    f"Sample: {sample_ids[j]}",
                    f"Population: {pop_labels[j]}",
                    f"PC1: {pcs_arr[j, 0]:.4f}",
                    f"PC2: {pcs_arr[j, 1]:.4f}",
                    f"PC3: {pcs_arr[j, 2]:.4f}",
                ]
                hover_texts.append("<br>".join(parts))

            fig.add_trace(
                go.Scatter3d(
                    x=pcs_arr[indices, 0].tolist(),
                    y=pcs_arr[indices, 1].tolist(),
                    z=pcs_arr[indices, 2].tolist(),
                    mode="markers",
                    marker=dict(
                        size=5,
                        color=palette[i % len(palette)],
                        opacity=0.8,
                    ),
                    text=hover_texts,
                    hoverinfo="text",
                    name=pop,
                )
            )

        fig.update_layout(
            title=title,
            scene=dict(
                xaxis_title=_axis_label(0),
                yaxis_title=_axis_label(1),
                zaxis_title=_axis_label(2),
            ),
            template="plotly_white",
        )

        if out is not None:
            fig.write_html(str(out), include_plotlyjs=True)
            logger.info(f"Saved interactive 3D PCA to {out}")

        return {
            "status": "success",
            "output_path": str(out) if out else None,
            "n_samples": n_samples,
            "interactive": True,
        }

    # --- Fallback: HTML summary table ----------------------------------
    logger.info("Plotly not available; generating fallback HTML table for PCA")
    headers = ["Sample", "Population", "PC1", "PC2", "PC3"]
    rows = []
    for i in range(n_samples):
        rows.append(
            [
                sample_ids[i],
                pop_labels[i],
                f"{pcs_arr[i, 0]:.4f}",
                f"{pcs_arr[i, 1]:.4f}",
                f"{pcs_arr[i, 2]:.4f}",
            ]
        )
    if out is not None:
        _fallback_html_table(title, headers, rows, out)
        logger.info(f"Saved fallback PCA table to {out}")

    return {
        "status": "success",
        "output_path": str(out) if out else None,
        "n_samples": n_samples,
        "interactive": False,
    }


# ---------------------------------------------------------------------------
# 3. Interactive Volcano Plot
# ---------------------------------------------------------------------------


def interactive_volcano(
    assoc_results: List[Dict[str, Any]],
    output_file: Optional[Union[str, Path]] = None,
    significance_threshold: float = 5e-8,
    fold_change_threshold: float = 0.5,
    title: str = "Interactive Volcano Plot",
) -> Dict[str, Any]:
    """Create an interactive volcano plot of GWAS association results.

    Generates a Plotly scatter with beta/effect_size on the x-axis and
    -log10(p_value) on the y-axis.  Points are coloured by significance
    category.  Falls back to a styled HTML table when Plotly is unavailable.

    Args:
        assoc_results: List of dicts, each with ``"p_value"`` (or ``"pvalue"``)
            and ``"beta"`` (or ``"effect_size"``).  Optional: ``"variant_id"``,
            and the chromosome (``"chromosome"``, ``"chrom"``, ``"chr"``,
            ``"CHROM"`` or ``"CHR"``).
        output_file: Path for the output HTML file.
        significance_threshold: P-value threshold for genome-wide significance.
        fold_change_threshold: Effect-size threshold for vertical guide lines.
        title: Plot title.

    Returns:
        Dict with ``"status"``, ``"output_path"``, ``"n_variants"``, and
        ``"interactive"`` keys.
    """
    if not assoc_results:
        logger.warning("No association results provided for volcano plot")
        return {
            "status": "skipped",
            "output_path": None,
            "n_variants": 0,
            "interactive": False,
        }

    out = _ensure_output_dir(output_file)

    # --- Normalise keys ------------------------------------------------
    records: List[Dict[str, Any]] = []
    for r in assoc_results:
        pval = r.get("p_value", r.get("pvalue", r.get("P", None)))
        beta = r.get("beta", r.get("effect_size", 0.0))
        if pval is None or pval <= 0:
            continue
        vid = r.get("variant_id", r.get("rsid", ""))
        chrom = _chromosome_label(_first_present(r, _CHROMOSOME_KEYS)) or ""
        neg_log_p = -math.log10(max(pval, 1e-300))

        # Categorise
        if pval < significance_threshold and abs(beta) >= fold_change_threshold:
            category = "Significant"
            colour = "red"
        elif pval < significance_threshold:
            category = "Suggestive"
            colour = "orange"
        else:
            category = "Not significant"
            colour = "gray"

        records.append(
            {
                "variant_id": vid,
                "chromosome": chrom,
                "p_value": float(pval),
                "beta": float(beta),
                "neg_log_p": neg_log_p,
                "category": category,
                "colour": colour,
            }
        )

    if not records:
        logger.warning("No valid variants after filtering for volcano plot")
        return {
            "status": "skipped",
            "output_path": None,
            "n_variants": 0,
            "interactive": False,
        }

    n_variants = len(records)

    # --- Plotly path ---------------------------------------------------
    if HAS_PLOTLY and go is not None:
        fig = go.Figure()

        category_colours = {
            "Significant": "red",
            "Suggestive": "orange",
            "Not significant": "gray",
        }

        for cat, col in category_colours.items():
            cat_records = [r for r in records if r["category"] == cat]
            if not cat_records:
                continue

            hover_texts = []
            for r in cat_records:
                parts = [
                    f"Variant: {r['variant_id']}" if r["variant_id"] else "",
                    f"Beta: {r['beta']:.4f}",
                    f"P-value: {r['p_value']:.2e}",
                ]
                hover_texts.append("<br>".join(p for p in parts if p))

            fig.add_trace(
                go.Scattergl(
                    x=[r["beta"] for r in cat_records],
                    y=[r["neg_log_p"] for r in cat_records],
                    mode="markers",
                    marker=dict(color=col, size=5, opacity=0.7),
                    text=hover_texts,
                    hoverinfo="text",
                    name=cat,
                )
            )

        # Threshold lines
        sig_y = -math.log10(significance_threshold)
        fig.add_hline(
            y=sig_y,
            line_dash="dash",
            line_color="red",
            annotation_text=f"p = {significance_threshold:.0e}",
        )
        fig.add_vline(x=fold_change_threshold, line_dash="dash", line_color="blue")
        fig.add_vline(x=-fold_change_threshold, line_dash="dash", line_color="blue")

        fig.update_layout(
            title=title,
            xaxis_title="Effect Size (Beta)",
            yaxis_title="-log<sub>10</sub>(p-value)",
            hovermode="closest",
            template="plotly_white",
        )

        if out is not None:
            fig.write_html(str(out), include_plotlyjs=True)
            logger.info(f"Saved interactive volcano plot to {out}")

        return {
            "status": "success",
            "output_path": str(out) if out else None,
            "n_variants": n_variants,
            "interactive": True,
        }

    # --- Fallback: styled HTML table -----------------------------------
    logger.info("Plotly not available; generating fallback HTML table for volcano plot")
    sorted_records = sorted(records, key=lambda r: r["p_value"])
    headers = ["Variant", "Chromosome", "Beta", "P-value", "-log10(p)", "Category"]
    rows = []
    for r in sorted_records:
        rows.append(
            [
                r["variant_id"] or "N/A",
                r["chromosome"] or "N/A",
                f"{r['beta']:.4f}",
                f"{r['p_value']:.2e}",
                f"{r['neg_log_p']:.2f}",
                r["category"],
            ]
        )
    if out is not None:
        _fallback_html_table(title, headers, rows, out)
        logger.info(f"Saved fallback volcano table to {out}")

    return {
        "status": "success",
        "output_path": str(out) if out else None,
        "n_variants": n_variants,
        "interactive": False,
    }
