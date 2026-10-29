"""Tests for metainformant.popgen.workflow.reporting (zero-mocks).

The report is generated from real seeded simulations analysed by the real scenario pipeline, written to a real
temporary directory and read back.
"""

from __future__ import annotations

import random
import re
from pathlib import Path
from typing import Any

import pytest

from metainformant.core.io import dump_json
from metainformant.popgen.workflow import analyze_dataset, generate_summary_report
from metainformant.simulation.models.popgen import (
    generate_genotype_matrix,
    generate_population_sequences,
    generate_two_populations,
    simulate_bottleneck_population,
    simulate_population_expansion,
)


def _write_fasta(sequences: list[str], path: Path) -> str:
    path.write_text("".join(f">seq_{i}\n{seq}\n" for i, seq in enumerate(sequences)), encoding="utf-8")
    return str(path)


@pytest.fixture(scope="module")
def report(tmp_path_factory: pytest.TempPathFactory) -> tuple[str, dict[str, Any]]:
    out = tmp_path_factory.mktemp("popgen_report")
    rng = random.Random(21)
    n_sequences, length = 14, 400
    scenarios: dict[str, dict[str, Any]] = {}
    for name, sequences in (
        ("neutral", generate_population_sequences(n_sequences, length, nucleotide_diversity=0.02, rng=rng)),
        (
            "bottleneck",
            simulate_bottleneck_population(
                n_sequences=n_sequences, sequence_length=length, bottleneck_size=5, bottleneck_duration=10, rng=rng
            ),
        ),
        (
            "expansion",
            simulate_population_expansion(
                n_sequences=n_sequences, sequence_length=length, expansion_factor=100.0, expansion_time=400, rng=rng
            ),
        ),
    ):
        scenarios[name] = {"file": _write_fasta(list(sequences), out / f"{name}.fasta")}
    for name, fst in (("two_populations_low_fst", 0.05), ("two_populations_high_fst", 0.3)):
        pop1, pop2 = generate_two_populations(10, 10, length, n_sites=length, fst=fst, rng=rng)
        scenarios[name] = {
            "file_pop1": _write_fasta(pop1, out / f"{name}_1.fasta"),
            "file_pop2": _write_fasta(pop2, out / f"{name}_2.fasta"),
        }
    genotypes = generate_genotype_matrix(30, 50, maf_min=0.1, maf_max=0.5, rng=rng)
    dump_json(genotypes, str(out / "genotypes.json"))
    scenarios["large_genotypes"] = {"file": str(out / "genotypes.json"), "n_individuals": 30, "n_sites": 50}
    dataset_info = {
        "generation_time": "2026-10-28T00:00:00+00:00",
        "seed": 21,
        "n_sequences_per_scenario": n_sequences,
        "sequence_length": length,
        "scenarios": scenarios,
    }
    results = analyze_dataset(dataset_info, out)
    generate_summary_report(dataset_info, results, out)
    return (out / "analysis_report.md").read_text(encoding="utf-8"), results


def test_the_conclusions_quote_the_tajima_d_of_each_scenario_from_the_analysis(report: tuple[str, Any]) -> None:
    text, results = report
    comparison = results["comparative_analysis"]["tajimas_d_comparison"]
    line = next(line for line in text.splitlines() if "Demographic signatures" in line)
    quoted = {name: float(value) for name, value in re.findall(r"(neutral|bottleneck|expansion) (-?\d+\.\d+)", line)}
    assert set(quoted) == {"neutral", "bottleneck", "expansion"}
    for name, value in quoted.items():
        assert value == pytest.approx(comparison[name], abs=5e-4)


def test_the_report_no_longer_claims_that_bottlenecks_and_expansions_both_give_negative_d(
    report: tuple[str, Any],
) -> None:
    text, _ = report
    assert "(negative Tajima's D)" not in text
    assert "Successfully generated populations with target diversity levels" not in text
    assert "Fst values match target specifications" not in text
    assert "A bottleneck sampled soon afterwards lacks rare variants and has a D above the neutral value" in text


def test_the_conclusions_quote_diversity_fst_and_the_genotype_matrix_of_the_dataset(
    report: tuple[str, Any],
) -> None:
    text, results = report
    comparative = results["comparative_analysis"]
    conclusions = text[text.index("## Conclusions") :]
    for name, value in comparative["diversity_comparison"].items():
        assert f"{name.replace('_', ' ')} {value:.6f}" in conclusions
    for name, value in comparative["fst_comparison"].items():
        assert f"{name.replace('_', ' ')} {value:.4f}" in conclusions
    assert "a genotype matrix of 30 individuals × 50 sites" in conclusions
