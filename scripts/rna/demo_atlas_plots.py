"""Render atlas demo figures from small synthetic data into docs figure set.

Thin orchestrator for the atlas_plots demo; descriptive-only.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from metainformant.rna.analysis.atlas_plots import (
    compute_tau,
    plot_orthogroup_small_multiples,
    plot_tau_heatmap,
    plot_tau_orthology_strips,
)

OUT = Path("docs/rna/figures/atlas_demo")
RNG = np.random.default_rng(20260901)

SPECIES = [
    "apis_mellifera",
    "bombus_terrestris",
    "ceratina_calcarata",
    "apis_dorsata",
    "melipona_quadrifasciata",
]
TISSUES = ["antenna", "brain", "midgut", "ovary", "testes"]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # (a) species x tissue mean-tau heatmap
    species_profiles = {sp: RNG.uniform(0.1, 0.95, size=len(TISSUES)) for sp in SPECIES}
    tau_frame = pd.DataFrame(
        {
            t: [species_profiles[sp][j] + RNG.normal(0, 0.03) for sp in SPECIES]
            for j, t in enumerate(TISSUES)
        },
        index=SPECIES,
    ).clip(0, 1)
    plot_tau_heatmap(tau_frame, OUT / "tau_heatmap.png", annot=True)

    # (b) per-orthogroup small multiples from synthetic long expression
    rows = []
    for og in [f"OG{i:04d}" for i in range(1, 7)]:
        base = float(RNG.uniform(10, 5000))
        for sp in SPECIES[:4]:
            for t in TISSUES:
                rows.append(
                    {
                        "orthogroup": og,
                        "species": sp,
                        "tissue": t,
                        "expression": max(base * RNG.uniform(0.05, 6.0), 1e-2),
                    }
                )
    plot_orthogroup_small_multiples(pd.DataFrame(rows), OUT / "orthogroup_profiles.png")

    # (c) tau vs orthology class strips, derived from synthetic per-gene tau
    class_shift = {
        "one_to_one": 0.25,
        "one_to_many": 0.15,
        "many_to_one": 0.05,
        "many_to_many": -0.05,
    }
    frames = []
    for cls, shift in class_shift.items():
        tau_vals = compute_tau(
            pd.DataFrame(
                RNG.dirichlet(np.full(len(TISSUES), 1.0 + shift + 0.4), size=150),
                index=[f"g_{cls}_{i}" for i in range(150)],
                columns=TISSUES,
            )
        )
        frames.append(pd.DataFrame({"tau": tau_vals.values, "orthology_class": cls}))
    plot_tau_orthology_strips(
        pd.concat(frames, ignore_index=True), OUT / "tau_by_orthology_class.png"
    )

    for f in sorted(OUT.glob("*.png")):
        print(f, f.stat().st_size)


if __name__ == "__main__":
    main()
