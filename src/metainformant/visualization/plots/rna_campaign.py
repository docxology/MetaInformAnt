from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

from metainformant.visualization.config.conventions import (
    CATEGORY_HATCHES,
    OKABE_ITO,
    save_figure_deterministic,
)


def plot_acquisition_observation(
    payload: dict[str, Any], output_prefix: Path, *, title: str = "RNA-seq acquisition"
) -> list[Path]:
    counts, scope = payload["counts"], payload["scope"]
    categories = [
        "Current lower bound",
        "Supported, unconfirmed",
        "Unresolved prerequisites",
    ]
    values = [
        counts["current_lower_bound"],
        counts["supported_unconfirmed"],
        scope["unsupported"],
    ]
    if any(type(value) is not int or value < 0 for value in values) or sum(values) != scope["total"]:
        raise ValueError("cohort chart partitions must equal the frozen total")
    phase_labels = [
        "Reused successes",
        "New, restore-bound",
        "New, unconfirmed",
        "Deferred downloads",
        "Failed, not deferred",
    ]
    phase_values = [
        counts["reused"],
        counts["new_restore_bound"],
        counts["new_unconfirmed"],
        counts["deferred_attempts"],
        counts["failed_nondeferred_attempts"],
    ]
    if any(type(value) is not int or value < 0 for value in phase_values):
        raise ValueError("phase counts must be nonnegative integers")
    if sum(phase_values[3:]) != counts["failed_attempts"]:
        raise ValueError("deferred and failed attempts must equal unsuccessful attempts")
    observed = datetime.fromtimestamp(payload["observed_at_epoch"], UTC).strftime("%Y-%m-%d %H:%M:%S UTC")
    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    ):
        figure, axes = plt.subplots(1, 2, figsize=(11.8, 4.8), gridspec_kw={"width_ratios": [1.35, 1]})
        for index, (label, value) in enumerate(zip(categories, values, strict=True)):
            bar = axes[0].barh(
                index,
                value,
                color=OKABE_ITO[(2, 4, 0)[index]],
                hatch=CATEGORY_HATCHES[index],
                edgecolor="white",
                height=0.58,
            )
            axes[0].bar_label(bar, labels=[f"{value:,}"], padding=6)
        axes[0].set_yticks(range(3), categories)
        axes[0].invert_yaxis()
        axes[0].set_xlim(0, max(values, default=1) * 1.18 + 1)
        axes[0].set_xlabel("Run accessions")
        axes[0].set_title(
            f"A  Frozen inventory: {scope['total']:,} runs",
            loc="left",
            fontweight="bold",
        )
        bars = axes[1].barh(
            range(5),
            phase_values,
            color=[OKABE_ITO[i] for i in (4, 2, 0, 3, 5)],
            height=0.58,
        )
        for index, bar in enumerate(bars):
            bar.set_hatch(CATEGORY_HATCHES[index])
        axes[1].bar_label(bars, labels=[f"{value:,}" for value in phase_values], padding=6)
        axes[1].set_yticks(range(5), phase_labels)
        axes[1].invert_yaxis()
        axes[1].set_xlim(0, max(phase_values, default=1) * 1.25 + 1)
        axes[1].set_xlabel("Terminal journal attempts")
        axes[1].set_title("B  Explicit phase journal", loc="left", fontweight="bold")
        figure.suptitle(f"{title} | operational observation", x=0.025, ha="left", fontsize=14)
        transfer = payload["transfer"]
        rate = (
            f"Retained-growth lower bound: {transfer['MiB_s_lower_bound']:.2f} MiB/s over {transfer['seconds']:.1f} s"
            if transfer["available"]
            else f"Transfer rate unavailable: {transfer['reason']}"
        )
        figure.text(0.025, 0.12, f"{observed}  |  {rate}", fontsize=9)
        figure.text(
            0.025,
            0.065,
            "Startup counts are producer-audited; new increments require matching restore evidence.",
            fontsize=9,
        )
        figure.text(
            0.025,
            0.025,
            "Operational counts are not biological replicates. This is not a scientific release.",
            fontsize=9,
        )
        figure.subplots_adjust(left=0.21, right=0.97, top=0.79, bottom=0.28, wspace=0.82)
        outputs = [output_prefix.with_suffix(suffix) for suffix in (".png", ".pdf")]
        try:
            for path in outputs:
                save_figure_deterministic(figure, path, dpi=180)
        finally:
            plt.close(figure)
    return outputs
