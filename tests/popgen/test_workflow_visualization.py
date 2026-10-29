"""Tests for metainformant.popgen.workflow.visualization.

Uses the real seeded simulation generators, the real scenario analysis and real
matplotlib rendering on a small dataset; nothing is replaced by a test double.
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any

import pytest

from metainformant.core.io import dump_json
from metainformant.popgen.workflow import analyze_dataset, generate_visualizations
from metainformant.simulation.models.popgen import (
    generate_genotype_matrix,
    generate_linkage_disequilibrium_data,
    generate_population_sequences,
    generate_two_populations,
    simulate_bottleneck_population,
    simulate_population_expansion,
)

_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"

# One plot per visualization step that has data in the dataset below.
_EXPECTED_PLOTS = {
    "diversity_comparison.png",
    "tajimas_d_comparison.png",
    "fst_comparison.png",
    "neutrality_test_summary.png",
    "pca_analysis.png",
    "kinship_matrix.png",
    "site_frequency_spectrum_example.png",
    "demographic_model_comparison.png",
    "summary_statistics_grid.png",
    "linkage_disequilibrium_decay.png",
    "comprehensive_neutrality_test_suite.png",
    "statistic_correlation_matrix.png",
}


def _write_fasta(sequences: list[str], path: Path) -> str:
    path.write_text("".join(f">seq_{i}\n{seq}\n" for i, seq in enumerate(sequences)), encoding="utf-8")
    return str(path)


def _small_dataset_info(data_dir: Path) -> dict[str, Any]:
    """Build the dataset description used by analyze_dataset from small real simulations."""
    rng = random.Random(11)
    n_sequences, length = 12, 240
    scenarios: dict[str, dict[str, Any]] = {}

    for name, diversity, mutation_rate in (
        ("neutral", 0.01, 0.001),
        ("high_diversity", 0.05, 0.005),
        ("low_diversity", 0.001, 0.0001),
    ):
        sequences = generate_population_sequences(
            n_sequences=n_sequences,
            sequence_length=length,
            nucleotide_diversity=diversity,
            mutation_rate=mutation_rate,
            rng=rng,
        )
        scenarios[name] = {"file": _write_fasta(sequences, data_dir / f"{name}.fasta")}

    bottleneck = simulate_bottleneck_population(
        n_sequences=n_sequences,
        sequence_length=length,
        pre_bottleneck_diversity=0.01,
        bottleneck_size=5,
        bottleneck_duration=10,
        mutation_rate=0.001,
        rng=rng,
    )
    scenarios["bottleneck"] = {"file": _write_fasta(bottleneck, data_dir / "bottleneck.fasta")}

    expansion = simulate_population_expansion(
        n_sequences=n_sequences,
        sequence_length=length,
        expansion_factor=10.0,
        mutation_rate=0.001,
        rng=rng,
    )
    scenarios["expansion"] = {"file": _write_fasta(expansion, data_dir / "expansion.fasta")}

    for name, fst in (("two_populations_low_fst", 0.05), ("two_populations_high_fst", 0.3)):
        pop1, pop2 = generate_two_populations(
            n_pop1=10, n_pop2=10, sequence_length=length, n_sites=length, fst=fst, rng=rng
        )
        scenarios[name] = {
            "file_pop1": _write_fasta(pop1, data_dir / f"{name}_pop1.fasta"),
            "file_pop2": _write_fasta(pop2, data_dir / f"{name}_pop2.fasta"),
        }

    genotypes = generate_genotype_matrix(n_individuals=40, n_sites=120, maf_min=0.05, maf_max=0.5, hwe=True, rng=rng)
    genotype_file = data_dir / "genotypes.json"
    dump_json(genotypes, str(genotype_file))
    scenarios["large_genotypes"] = {"file": str(genotype_file)}

    ld_genotypes = generate_linkage_disequilibrium_data(n_individuals=40, n_sites=60, r_squared_target=0.5, rng=rng)
    ld_file = data_dir / "ld_genotypes.json"
    dump_json(ld_genotypes, str(ld_file))
    scenarios["linkage_disequilibrium"] = {"file": str(ld_file)}

    return {"scenarios": scenarios}


@pytest.fixture(scope="module")
def plots_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    output_dir = tmp_path_factory.mktemp("popgen_workflow")
    data_dir = output_dir / "data"
    data_dir.mkdir()
    analysis_results = analyze_dataset(_small_dataset_info(data_dir), output_dir)
    assert analysis_results["validation"]["status"] == "passed"
    assert generate_visualizations(analysis_results, output_dir) is None
    return output_dir / "plots"


class TestGenerateVisualizations:
    def test_every_plot_step_writes_a_png(self, plots_dir: Path) -> None:
        produced = {path.name for path in plots_dir.iterdir()}
        assert _EXPECTED_PLOTS <= produced

    @pytest.mark.parametrize("name", sorted(_EXPECTED_PLOTS))
    def test_plot_is_a_nonempty_png_file(self, plots_dir: Path, name: str) -> None:
        content = (plots_dir / name).read_bytes()
        assert content.startswith(_PNG_SIGNATURE)
        assert len(content) > 1000

    def test_summary_fst_values_do_not_produce_an_fst_matrix(self, plots_dir: Path) -> None:
        # plot_fst_matrix needs population->sequences input; the scenario analyses only hold
        # summary Fst values, so the workflow deliberately skips that plot.
        assert not (plots_dir / "fst_matrix.png").exists()
        assert {path.name for path in plots_dir.iterdir()} == _EXPECTED_PLOTS
