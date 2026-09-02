"""Tests that scripts/popgen imports resolve against the new src module.

Guards the thin-orchestrator contract: every metainformant.* import used
by the scripts must exist (this regressed when
``metainformant.dna.population_analysis`` and
``metainformant.simulation.popgen`` were refactored away).
"""

from __future__ import annotations

import importlib

import pytest

SCRIPT_IMPORTS = [
    "metainformant.dna.population.analysis",  # calculate_summary_statistics etc.
    "metainformant.dna.population.core",  # fu_li/fay_wu helpers
    "metainformant.dna.sequence.core",  # read_fasta
    "metainformant.gwas.analysis.quality",  # test_hwe
    "metainformant.gwas.analysis.structure",  # compute_pca / compute_kinship_matrix
    "metainformant.math.population_genetics.demography",
    "metainformant.math.population_genetics.ld",
    "metainformant.simulation.models.popgen",
    "metainformant.dna.population.visualization_stats",
    "metainformant.dna.population.visualization_core",
]


@pytest.mark.parametrize("module_name", SCRIPT_IMPORTS)
def test_script_dependency_module_imports(module_name: str) -> None:
    importlib.import_module(module_name)


def test_renamed_symbols_exist() -> None:
    """Symbols scripts/popgen previously imported from dead paths."""
    from metainformant.dna.population.analysis import (
        calculate_summary_statistics,
        compare_populations,
        neutrality_test_suite,
    )
    from metainformant.simulation.models.popgen import (
        generate_site_frequency_spectrum,
        generate_genotype_matrix,
        generate_population_sequences,
    )

    assert callable(calculate_summary_statistics)
    assert callable(compare_populations)
    assert callable(neutrality_test_suite)
    assert callable(generate_population_sequences)
    assert callable(generate_genotype_matrix)
    assert callable(generate_site_frequency_spectrum)


def test_popgen_module_exports() -> None:
    import metainformant.popgen as pg

    for name in (
        "analyze_dataset",
        "summarize_scenario",
        "sequence_scenario_suite",
        "compare_two_population_sequences",
        "genotype_structure_analysis",
        "ld_summary",
        "demographic_model_comparisons",
    ):
        assert hasattr(pg, name), f"missing export: {name}"
