"""Export the offline reader and regenerate its figure from a synthetic TSV."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import matplotlib
from matplotlib import pyplot as plt
from metainformant.core import io

PROFILE = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path, help="New reader project directory")
    destination = parser.parse_args().destination
    if destination.exists():
        parser.error("destination must be a new directory")

    rows = list(io.read_delimited(PROFILE / "synthetic_counts.tsv", delimiter="\t"))
    if not rows or any(list(row) != ["group", "count"] for row in rows):
        raise ValueError("Expected ordered group and count columns")
    groups = [row["group"] for row in rows]
    counts = [int(row["count"]) for row in rows]
    if len(set(groups)) != len(groups) or any(not group for group in groups):
        raise ValueError("Expected nonempty, unique group labels")
    if any(count < 0 for count in counts):
        raise ValueError("Expected nonnegative counts")

    source = destination / "docs/manuscript"
    data = destination / "output/data"
    figures = destination / "output/figures"
    for directory in (source, data, figures):
        directory.mkdir(parents=True)
    for name in ("01_scope.md", "02_table_check.md", "config.yaml", "references.bib"):
        shutil.copy2(PROFILE / name, source / name)

    roundtrip = data / "synthetic_counts.tsv"
    io.write_delimited(rows, roundtrip, delimiter="\t")
    if list(io.read_delimited(roundtrip, delimiter="\t")) != rows:
        raise ValueError("Synthetic counts changed on round trip")

    (data / "synthetic_count_summary.json").write_text(
        json.dumps(dict(zip(groups, counts, strict=True)), indent=2) + "\n",
        encoding="utf-8",
    )

    matplotlib.use("Agg")
    fig, ax = plt.subplots(figsize=(6, 3.5), layout="constrained")
    bars = ax.bar(groups, counts, color="#27627a")
    ax.bar_label(bars, padding=3)
    ax.set(xlabel="Invented group", ylabel="Count", ylim=(0, max(counts) + 1))
    fig.savefig(figures / "synthetic_table_counts.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
