"""Pipeline progress dashboard — mosaic graphical abstract.

Generates a comprehensive PDF + PNG dashboard showing sample processing status
across all species. Designed for quick visual assessment of pipeline health.

Usage::

    from metainformant.rna.engine.progress_dashboard import generate_dashboard
    generate_dashboard()  # writes to output/amalgkit/pipeline_dashboard.{pdf,png}

Or via CLI::

    uv run python scripts/rna/check_pipeline_status.py --dashboard
"""

from __future__ import annotations

import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import numpy as np

# ═══════════════════════════════════════════════════════════════
# Color palette — clean, print-friendly on white
# ═══════════════════════════════════════════════════════════════

COLORS = {
    "pending":     "#94a3b8",  # cool gray
    "downloading": "#3b82f6",  # vivid blue
    "downloaded":  "#06b6d4",  # cyan
    "quantifying": "#f59e0b",  # amber
    "quantified":  "#10b981",  # emerald
    "failed":      "#ef4444",  # red
}

STATE_ORDER = ["pending", "downloading", "downloaded", "quantifying", "quantified", "failed"]
STATE_LABELS = {
    "pending": "Pending",
    "downloading": "Downloading",
    "downloaded": "Downloaded",
    "quantifying": "Quantifying",
    "quantified": "Quantified",
    "failed": "Failed",
}

BG_COLOR = "#ffffff"
CARD_BG = "#f8fafc"
TEXT_COLOR = "#1e293b"
ACCENT = "#0f172a"
ACCENT2 = "#0369a1"       # secondary accent (blue)
GRID_COLOR = "#e2e8f0"
BORDER_COLOR = "#cbd5e1"
MUTED = "#64748b"

SPECIES_ORDER = [
    "anoplolepis_gracilipes", "acromyrmex_echinatior", "dinoponera_quadriceps",
    "vollenhovia_emeryi", "odontomachus_brunneus", "formica_exsecta",
    "temnothorax_americanus", "wasmannia_auropunctata", "nylanderia_fulva",
    "temnothorax_curvispinosus", "pbarbatus", "cardiocondyla_obscurior",
    "temnothorax_nylanderi", "linepithema_humile", "atta_cephalotes",
    "ooceraea_biroi", "camponotus_floridanus", "solenopsis_invicta",
    "monomorium_pharaonis", "temnothorax_longispinosus", "harpegnathos_saltator",
    "amellifera",
]

# Nicer display names (italic genus abbreviation + species)
def _display_name(sp: str) -> str:
    parts = sp.replace("_", " ").split()
    if len(parts) >= 2:
        return f"{parts[0][0].upper()}. {' '.join(p.capitalize() for p in parts[1:])}"
    return sp.replace("_", " ").title()

DEFAULT_DB = Path("output/amalgkit/pipeline_progress.db")
DEFAULT_OUTPUT = Path("output/amalgkit")


# ═══════════════════════════════════════════════════════════════
# Data loading
# ═══════════════════════════════════════════════════════════════

def load_counts(db_path: Path = DEFAULT_DB) -> Dict[str, Dict[str, int]]:
    """Load species × state counts from the progress database."""
    conn = sqlite3.connect(str(db_path))
    rows = conn.execute(
        "SELECT species, state, COUNT(*) FROM samples GROUP BY species, state"
    ).fetchall()
    conn.close()

    result: Dict[str, Dict[str, int]] = {}
    for sp, state, count in rows:
        if sp not in result:
            result[sp] = {}
        result[sp][state] = count
    return result


def load_failed_details(db_path: Path = DEFAULT_DB) -> List[Dict[str, Any]]:
    """Load failed sample details."""
    conn = sqlite3.connect(str(db_path))
    rows = conn.execute(
        "SELECT species, srr_id, error, updated_at FROM samples "
        "WHERE state = 'failed' ORDER BY species, updated_at DESC"
    ).fetchall()
    conn.close()
    return [{"species": r[0], "srr_id": r[1], "error": r[2], "when": r[3]} for r in rows]


# ═══════════════════════════════════════════════════════════════
# Individual plot components
# ═══════════════════════════════════════════════════════════════

def _style_ax(ax: plt.Axes, title: str = "") -> None:
    """Apply consistent styling to an axes."""
    ax.set_facecolor(CARD_BG)
    for spine in ax.spines.values():
        spine.set_color(BORDER_COLOR)
        spine.set_linewidth(0.8)
    ax.tick_params(colors=TEXT_COLOR, labelsize=10)
    if title:
        ax.set_title(title, color=ACCENT, fontsize=14, fontweight="bold",
                      pad=12, loc="left")


def plot_species_bars(ax: plt.Axes, counts: Dict[str, Dict[str, int]]) -> None:
    """Horizontal stacked bar chart — all species, all states."""
    species = [sp for sp in SPECIES_ORDER if sp in counts]
    if not species:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", color=MUTED, fontsize=14)
        return

    bar_data = {s: [] for s in STATE_ORDER}
    for sp in species:
        for state in STATE_ORDER:
            bar_data[state].append(counts.get(sp, {}).get(state, 0))

    y_pos = list(range(len(species)))
    left = [0] * len(species)

    for state in STATE_ORDER:
        vals = bar_data[state]
        if any(v > 0 for v in vals):
            ax.barh(y_pos, vals, left=left, height=0.65,
                    color=COLORS[state], label=STATE_LABELS[state],
                    edgecolor="white", linewidth=0.3)
            left = [l + v for l, v in zip(left, vals)]

    # Total count at end of each bar
    for i, sp in enumerate(species):
        total = sum(counts.get(sp, {}).values())
        ax.text(left[i] + max(left) * 0.01, i, f" {total:,}",
                va="center", ha="left", color=MUTED, fontsize=8, fontweight="bold")

    display_names = [_display_name(sp) for sp in species]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(display_names, fontsize=9, color=TEXT_COLOR)
    ax.invert_yaxis()
    ax.set_xlabel("Number of Samples", color=TEXT_COLOR, fontsize=11)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
    _style_ax(ax, "Samples per Species by Processing State")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_completion_pct(ax: plt.Axes, counts: Dict[str, Dict[str, int]]) -> None:
    """Horizontal completion percentage bars with color coding."""
    species = [sp for sp in SPECIES_ORDER if sp in counts]
    if not species:
        return

    pcts = []
    bg_colors = []
    fg_colors = []
    for sp in species:
        sc = counts.get(sp, {})
        total = sum(sc.values())
        quant = sc.get("quantified", 0)
        failed = sc.get("failed", 0)
        pct = (quant / total * 100) if total > 0 else 0
        pcts.append(pct)

        if pct >= 80:
            bg_colors.append("#10b981")
            fg_colors.append("#065f46")
        elif pct >= 20:
            bg_colors.append("#f59e0b")
            fg_colors.append("#92400e")
        elif pct > 0:
            bg_colors.append("#3b82f6")
            fg_colors.append("#1e40af")
        else:
            bg_colors.append("#e2e8f0")
            fg_colors.append("#94a3b8")

    y_pos = list(range(len(species)))

    # Background track
    ax.barh(y_pos, [100] * len(species), height=0.65,
            color="#f1f5f9", edgecolor="none")

    # Foreground progress
    bars = ax.barh(y_pos, pcts, height=0.65,
                   color=bg_colors, edgecolor="white", linewidth=0.3)

    # Percentage labels
    for i, pct in enumerate(pcts):
        label_x = max(pct + 2, 8)
        ax.text(label_x, i, f"{pct:.0f}%",
                va="center", ha="left", color=fg_colors[i],
                fontsize=9, fontweight="bold")

    ax.set_yticks([])
    ax.set_xlim(0, 110)
    ax.set_xlabel("% Quantified", color=TEXT_COLOR, fontsize=11)
    ax.invert_yaxis()
    _style_ax(ax, "Completion")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


def plot_overall_donut(ax: plt.Axes, counts: Dict[str, Dict[str, int]]) -> None:
    """Overall pipeline donut chart with center stats."""
    totals = {s: 0 for s in STATE_ORDER}
    for sp_counts in counts.values():
        for state, count in sp_counts.items():
            if state in totals:
                totals[state] += count

    grand_total = sum(totals.values())
    if grand_total == 0:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", color=MUTED, fontsize=14)
        return

    labels = []
    sizes = []
    clrs = []
    for s in STATE_ORDER:
        if totals[s] > 0:
            labels.append(STATE_LABELS[s])
            sizes.append(totals[s])
            clrs.append(COLORS[s])

    wedges, _ = ax.pie(
        sizes, colors=clrs, startangle=90,
        wedgeprops={"width": 0.32, "edgecolor": "white", "linewidth": 2.5}
    )

    # Center stats
    quant = totals.get("quantified", 0)
    pct = quant / grand_total * 100
    ax.text(0, 0.08, f"{pct:.1f}%", ha="center", va="center",
            color=ACCENT, fontsize=28, fontweight="bold")
    ax.text(0, -0.12, f"{quant:,} / {grand_total:,}", ha="center", va="center",
            color=MUTED, fontsize=11)
    ax.text(0, -0.25, "quantified", ha="center", va="center",
            color=MUTED, fontsize=9, style="italic")

    ax.set_facecolor(CARD_BG)
    _style_ax(ax, "Overall Progress")


def plot_state_summary_table(ax: plt.Axes, counts: Dict[str, Dict[str, int]]) -> None:
    """Clean table showing aggregate state counts."""
    totals = {s: 0 for s in STATE_ORDER}
    for sp_counts in counts.values():
        for state, count in sp_counts.items():
            if state in totals:
                totals[state] += count

    grand_total = sum(totals.values())
    n_species = len(counts)
    n_complete = sum(
        1 for sp in counts
        if counts[sp].get("quantified", 0) == sum(counts[sp].values()) and sum(counts[sp].values()) > 0
    )

    ax.axis("off")
    ax.set_facecolor(CARD_BG)
    _style_ax(ax, "Pipeline Summary")

    y = 0.82
    line_h = 0.085

    def _row(label: str, value: str, color: str = TEXT_COLOR, bold: bool = False):
        nonlocal y
        weight = "bold" if bold else "normal"
        ax.text(0.05, y, label, transform=ax.transAxes, fontsize=11,
                color=MUTED, va="top", fontfamily="sans-serif")
        ax.text(0.95, y, value, transform=ax.transAxes, fontsize=11,
                color=color, va="top", ha="right", fontweight=weight, fontfamily="monospace")
        y -= line_h

    _row("Species tracked", str(n_species), ACCENT, True)
    _row("Total samples", f"{grand_total:,}", ACCENT, True)
    _row("Species complete", f"{n_complete} / {n_species}", COLORS["quantified"])
    _row("", "")  # spacer

    for state in STATE_ORDER:
        count = totals[state]
        pct = f"({count/grand_total*100:.1f}%)" if grand_total > 0 else ""
        _row(STATE_LABELS[state], f"{count:>7,}  {pct}", COLORS[state])


def plot_species_heatmap(ax: plt.Axes, counts: Dict[str, Dict[str, int]]) -> None:
    """Heatmap grid: species × state, cell intensity by count."""
    species = [sp for sp in SPECIES_ORDER if sp in counts]
    if not species:
        return

    data = []
    for sp in species:
        row = [counts.get(sp, {}).get(s, 0) for s in STATE_ORDER]
        data.append(row)

    data_arr = np.array(data, dtype=float)
    # Log-scale the data for visibility
    log_data = np.log1p(data_arr)
    max_val = log_data.max() if log_data.max() > 0 else 1

    ax.set_facecolor(CARD_BG)

    n_sp = len(species)
    n_st = len(STATE_ORDER)

    cell_w = 1.0
    cell_h = 0.8
    pad = 0.08

    for i, sp in enumerate(species):
        for j, state in enumerate(STATE_ORDER):
            val = data_arr[i][j]
            intensity = log_data[i][j] / max_val if max_val > 0 else 0

            # Cell color: blend from card bg to state color
            base_rgb = matplotlib.colors.to_rgb(COLORS[state])
            bg_rgb = matplotlib.colors.to_rgb(CARD_BG)
            blended = tuple(bg + (c - bg) * intensity for bg, c in zip(bg_rgb, base_rgb))

            rect = mpatches.FancyBboxPatch(
                (j * (cell_w + pad), (n_sp - 1 - i) * (cell_h + pad)),
                cell_w, cell_h,
                boxstyle="round,pad=0.03",
                facecolor=blended,
                edgecolor="white", linewidth=0.5,
            )
            ax.add_patch(rect)

            # Cell text
            if val > 0:
                text_color = "white" if intensity > 0.5 else TEXT_COLOR
                ax.text(
                    j * (cell_w + pad) + cell_w / 2,
                    (n_sp - 1 - i) * (cell_h + pad) + cell_h / 2,
                    f"{int(val):,}",
                    ha="center", va="center", fontsize=7,
                    color=text_color, fontweight="bold",
                )

    # Axis labels
    ax.set_xlim(-pad, n_st * (cell_w + pad))
    ax.set_ylim(-pad, n_sp * (cell_h + pad))
    ax.set_xticks([j * (cell_w + pad) + cell_w / 2 for j in range(n_st)])
    ax.set_xticklabels([STATE_LABELS[s] for s in STATE_ORDER],
                       fontsize=9, color=TEXT_COLOR, rotation=30, ha="right")
    ax.set_yticks([(n_sp - 1 - i) * (cell_h + pad) + cell_h / 2 for i in range(n_sp)])
    ax.set_yticklabels([_display_name(sp) for sp in species], fontsize=8, color=TEXT_COLOR)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    _style_ax(ax, "Status Heatmap (log-scaled)")


def plot_failed_summary(ax: plt.Axes, counts: Dict[str, Dict[str, int]],
                        failed_details: List[Dict[str, Any]]) -> None:
    """Summary of failed samples by species."""
    ax.axis("off")
    ax.set_facecolor(CARD_BG)
    _style_ax(ax, "Failed Samples")

    # Aggregate failures by species
    failures_by_sp: Dict[str, int] = {}
    for sp_counts in counts.items():
        sp, sc = sp_counts
        failed = sc.get("failed", 0)
        if failed > 0:
            failures_by_sp[sp] = failed

    if not failures_by_sp:
        ax.text(0.5, 0.5, "No failures 🎉", transform=ax.transAxes,
                ha="center", va="center", color=COLORS["quantified"],
                fontsize=14, fontweight="bold")
        return

    y = 0.82
    total_failed = sum(failures_by_sp.values())
    ax.text(0.05, y, f"Total failed: {total_failed:,}",
            transform=ax.transAxes, fontsize=12, color=COLORS["failed"],
            fontweight="bold")
    y -= 0.1

    sorted_failures = sorted(failures_by_sp.items(), key=lambda x: x[1], reverse=True)
    for sp, count in sorted_failures[:8]:
        ax.text(0.08, y, f"• {_display_name(sp)}: {count:,}",
                transform=ax.transAxes, fontsize=10, color=TEXT_COLOR)
        y -= 0.08

    if len(sorted_failures) > 8:
        remaining = sum(c for _, c in sorted_failures[8:])
        ax.text(0.08, y, f"  + {len(sorted_failures) - 8} more ({remaining:,})",
                transform=ax.transAxes, fontsize=9, color=MUTED, style="italic")


def plot_sample_counts_table(ax: plt.Axes, counts: Dict[str, Dict[str, int]]) -> None:
    """Detailed table with exact counts per species."""
    species = [sp for sp in SPECIES_ORDER if sp in counts]
    if not species:
        return

    ax.axis("off")
    ax.set_facecolor(CARD_BG)
    _style_ax(ax, "Detailed Sample Counts")

    # Headers
    headers = ["Species", "Total", "Quant", "Fail", "%"]
    col_x = [0.02, 0.45, 0.60, 0.75, 0.88]
    y = 0.90
    for h, x in zip(headers, col_x):
        ax.text(x, y, h, transform=ax.transAxes, fontsize=9,
                color=ACCENT, fontweight="bold", va="top")

    y -= 0.03
    ax.plot([0.02, 0.98], [y, y], color=BORDER_COLOR,
            linewidth=0.8, transform=ax.transAxes, clip_on=False)
    y -= 0.015

    line_h = 0.038 if len(species) > 15 else 0.042
    fsize = 7.5 if len(species) > 15 else 8.5

    for sp in species:
        sc = counts.get(sp, {})
        total = sum(sc.values())
        quant = sc.get("quantified", 0)
        failed = sc.get("failed", 0)
        pct = f"{quant/total*100:.0f}%" if total > 0 else "—"

        name_color = COLORS["quantified"] if quant == total and total > 0 else TEXT_COLOR

        ax.text(col_x[0], y, _display_name(sp), transform=ax.transAxes,
                fontsize=fsize, color=name_color, va="top")
        ax.text(col_x[1], y, f"{total:,}", transform=ax.transAxes,
                fontsize=fsize, color=TEXT_COLOR, va="top", fontfamily="monospace")
        ax.text(col_x[2], y, f"{quant:,}", transform=ax.transAxes,
                fontsize=fsize, color=COLORS["quantified"] if quant > 0 else MUTED, va="top",
                fontfamily="monospace")
        ax.text(col_x[3], y, f"{failed:,}", transform=ax.transAxes,
                fontsize=fsize, color=COLORS["failed"] if failed > 0 else MUTED, va="top",
                fontfamily="monospace")
        ax.text(col_x[4], y, pct, transform=ax.transAxes,
                fontsize=fsize, color=COLORS["quantified"] if quant > 0 else MUTED, va="top",
                fontweight="bold" if quant == total and total > 0 else "normal",
                fontfamily="monospace")
        y -= line_h

    # Totals row
    y -= 0.005
    ax.plot([0.02, 0.98], [y, y], color=BORDER_COLOR,
            linewidth=0.8, transform=ax.transAxes, clip_on=False)
    y -= 0.015
    grand_total = sum(sum(counts.get(sp, {}).values()) for sp in species)
    grand_quant = sum(counts.get(sp, {}).get("quantified", 0) for sp in species)
    grand_fail = sum(counts.get(sp, {}).get("failed", 0) for sp in species)
    grand_pct = f"{grand_quant/grand_total*100:.1f}%" if grand_total > 0 else "—"

    ax.text(col_x[0], y, "TOTAL", transform=ax.transAxes,
            fontsize=fsize + 1, color=ACCENT, va="top", fontweight="bold")
    ax.text(col_x[1], y, f"{grand_total:,}", transform=ax.transAxes,
            fontsize=fsize + 1, color=ACCENT, va="top", fontweight="bold", fontfamily="monospace")
    ax.text(col_x[2], y, f"{grand_quant:,}", transform=ax.transAxes,
            fontsize=fsize + 1, color=COLORS["quantified"], va="top", fontweight="bold",
            fontfamily="monospace")
    ax.text(col_x[3], y, f"{grand_fail:,}", transform=ax.transAxes,
            fontsize=fsize + 1, color=COLORS["failed"], va="top", fontweight="bold",
            fontfamily="monospace")
    ax.text(col_x[4], y, grand_pct, transform=ax.transAxes,
            fontsize=fsize + 1, color=ACCENT, va="top", fontweight="bold",
            fontfamily="monospace")


# ═══════════════════════════════════════════════════════════════
# Main dashboard composer
# ═══════════════════════════════════════════════════════════════

def generate_dashboard(
    db_path: Path = DEFAULT_DB,
    output_dir: Path = DEFAULT_OUTPUT,
    filename: str = "pipeline_dashboard",
) -> Tuple[Path, Path]:
    """Generate a comprehensive pipeline status dashboard.

    Args:
        db_path: Path to the SQLite progress database.
        output_dir: Directory to write PDF and PNG outputs.
        filename: Base filename (without extension).

    Returns:
        Tuple of (pdf_path, png_path).
    """
    counts = load_counts(db_path)
    failed_details = load_failed_details(db_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Figure: large, comprehensive
    fig = plt.figure(figsize=(24, 18), facecolor=BG_COLOR, dpi=150)

    # ── Title ──
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    grand_total = sum(sum(v.values()) for v in counts.values())
    grand_quant = sum(v.get("quantified", 0) for v in counts.values())

    fig.suptitle(
        "MetaInformAnt  ·  Amalgkit RNA-seq Pipeline Dashboard",
        fontsize=22, fontweight="bold", color=ACCENT,
        y=0.975, fontfamily="sans-serif"
    )
    fig.text(
        0.5, 0.96,
        f"Generated {timestamp}  ·  {len(counts)} species  ·  "
        f"{grand_total:,} total samples  ·  {grand_quant:,} quantified ({grand_quant/grand_total*100:.1f}%)"
        if grand_total > 0 else f"Generated {timestamp}  ·  No data",
        ha="center", fontsize=12, color=MUTED, fontfamily="sans-serif"
    )

    # ── Layout: 3 rows × 4 columns ──
    gs = gridspec.GridSpec(
        3, 4, figure=fig,
        left=0.03, right=0.98, top=0.94, bottom=0.03,
        hspace=0.30, wspace=0.22,
        height_ratios=[1.3, 1.0, 1.0],
    )

    # Row 1: Species bars (3 cols) + Completion % (1 col)
    ax_bars = fig.add_subplot(gs[0, :3])
    plot_species_bars(ax_bars, counts)

    ax_pct = fig.add_subplot(gs[0, 3])
    plot_completion_pct(ax_pct, counts)

    # Row 2: Donut (1 col) + Summary table (1 col) + Heatmap (2 cols)
    ax_donut = fig.add_subplot(gs[1, 0])
    plot_overall_donut(ax_donut, counts)

    ax_summary = fig.add_subplot(gs[1, 1])
    plot_state_summary_table(ax_summary, counts)

    ax_heatmap = fig.add_subplot(gs[1, 2:])
    plot_species_heatmap(ax_heatmap, counts)

    # Row 3: Detailed counts table (2 cols) + Failed summary (1 col) + Legend (1 col)
    ax_table = fig.add_subplot(gs[2, :2])
    plot_sample_counts_table(ax_table, counts)

    ax_failed = fig.add_subplot(gs[2, 2])
    plot_failed_summary(ax_failed, counts, failed_details)

    # Legend panel
    ax_legend = fig.add_subplot(gs[2, 3])
    ax_legend.set_facecolor(CARD_BG)
    ax_legend.axis("off")
    _style_ax(ax_legend, "Legend")

    legend_patches = [
        mpatches.Patch(color=COLORS[s], label=f"  {STATE_LABELS[s]}")
        for s in STATE_ORDER
    ]
    ax_legend.legend(
        handles=legend_patches, loc="center", ncol=1,
        fontsize=11, frameon=False, labelcolor=TEXT_COLOR,
        handlelength=2.5, handleheight=1.8,
    )

    # State flow diagram as text
    ax_legend.text(
        0.5, 0.12,
        "pending → downloading → downloaded\n"
        "→ quantifying → quantified | failed",
        transform=ax_legend.transAxes, fontsize=9,
        color=MUTED, ha="center", va="center",
        fontfamily="monospace", style="italic",
    )

    # ── Save ──
    pdf_path = output_dir / f"{filename}.pdf"
    png_path = output_dir / f"{filename}.png"

    fig.savefig(str(pdf_path), format="pdf", facecolor=BG_COLOR, bbox_inches="tight")
    fig.savefig(str(png_path), format="png", facecolor=BG_COLOR, bbox_inches="tight", dpi=200)
    plt.close(fig)

    print(f"  📊 Dashboard saved:")
    print(f"     PDF: {pdf_path}")
    print(f"     PNG: {png_path}")

    return pdf_path, png_path


if __name__ == "__main__":
    generate_dashboard()
