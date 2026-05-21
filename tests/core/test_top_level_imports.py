"""Top-level import and compatibility facade checks."""

from __future__ import annotations

import importlib

import pytest


TOP_LEVEL_PACKAGES = [
    "metainformant.cloud",
    "metainformant.core",
    "metainformant.dna",
    "metainformant.ecology",
    "metainformant.epigenome",
    "metainformant.gwas",
    "metainformant.information",
    "metainformant.life_events",
    "metainformant.longread",
    "metainformant.math",
    "metainformant.mcp",
    "metainformant.menu",
    "metainformant.metabolomics",
    "metainformant.metagenomics",
    "metainformant.ml",
    "metainformant.multiomics",
    "metainformant.networks",
    "metainformant.ontology",
    "metainformant.pharmacogenomics",
    "metainformant.phenotype",
    "metainformant.protein",
    "metainformant.quality",
    "metainformant.rna",
    "metainformant.simulation",
    "metainformant.singlecell",
    "metainformant.spatial",
    "metainformant.structural_variants",
    "metainformant.visualization",
]


@pytest.mark.parametrize("module_name", TOP_LEVEL_PACKAGES)
def test_top_level_package_imports(module_name: str) -> None:
    """Every top-level package advertised by README/CLI should import."""
    module = importlib.import_module(module_name)

    assert module.__name__ == module_name


def test_split_module_facades_preserve_public_imports() -> None:
    """Oversized modules keep their public import paths as facades."""
    from metainformant.gwas.visualization import _general_impl, general
    from metainformant.networks.interaction import _ppi_impl, ppi
    from metainformant.rna.amalgkit import _amalgkit_impl, amalgkit

    assert general.manhattan_plot is _general_impl.manhattan_plot
    assert amalgkit.run_amalgkit is _amalgkit_impl.run_amalgkit
    assert ppi.load_ppi_network is _ppi_impl.load_ppi_network


def test_gwas_top_level_convenience_export_is_lazy_accessible() -> None:
    """The lazy GWAS barrel preserves the historical convenience API."""
    from metainformant import gwas

    assert callable(gwas.association_test_linear)
