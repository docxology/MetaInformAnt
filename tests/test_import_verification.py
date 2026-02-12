"""Import verification tests for all METAINFORMANT modules.

Parametrized test that imports all 25 modules, verifies __all__ exports
are accessible, and checks canonical import paths.

NO MOCKING POLICY: All tests use real imports.
"""

from __future__ import annotations

import importlib

import pytest

# All 25 top-level modules
ALL_MODULES = [
    "metainformant.core",
    "metainformant.dna",
    "metainformant.rna",
    "metainformant.protein",
    "metainformant.gwas",
    "metainformant.math",
    "metainformant.information",
    "metainformant.life_events",
    "metainformant.visualization",
    "metainformant.networks",
    "metainformant.multiomics",
    "metainformant.singlecell",
    "metainformant.simulation",
    "metainformant.quality",
    "metainformant.ml",
    "metainformant.ontology",
    "metainformant.phenotype",
    "metainformant.ecology",
    "metainformant.epigenome",
    "metainformant.metagenomics",
    "metainformant.pharmacogenomics",
    "metainformant.spatial",
    "metainformant.structural_variants",
    "metainformant.longread",
    "metainformant.menu",
]


class TestModuleImports:
    """Verify all modules can be imported without errors."""

    @pytest.mark.parametrize("module_name", ALL_MODULES)
    def test_module_imports(self, module_name: str) -> None:
        """Each module should import without raising."""
        mod = importlib.import_module(module_name)
        assert mod is not None

    @pytest.mark.parametrize("module_name", ALL_MODULES)
    def test_module_has_all(self, module_name: str) -> None:
        """Each module should define __all__ for explicit exports."""
        mod = importlib.import_module(module_name)
        assert hasattr(mod, "__all__"), f"{module_name} is missing __all__"
        assert isinstance(mod.__all__, (list, tuple))
        assert len(mod.__all__) > 0, f"{module_name}.__all__ is empty"


class TestExportsAccessible:
    """Verify all __all__ exports are actually accessible."""

    @pytest.mark.parametrize("module_name", ALL_MODULES)
    def test_all_exports_accessible(self, module_name: str) -> None:
        """Every name in __all__ should be getattr-able from the module."""
        mod = importlib.import_module(module_name)
        if not hasattr(mod, "__all__"):
            pytest.skip(f"{module_name} has no __all__")

        missing = []
        for name in mod.__all__:
            if not hasattr(mod, name):
                missing.append(name)

        assert not missing, (
            f"{module_name} has {len(missing)} exports in __all__ that are not accessible: " f"{missing[:10]}"
        )


class TestCanonicalImportPaths:
    """Verify canonical import paths work for reorganized modules."""

    def test_ontology_canonical_paths(self) -> None:
        """Ontology canonical subpackage paths should work."""
        from metainformant.ontology.core import go, obo  # noqa: F401
        from metainformant.ontology.query import query  # noqa: F401
        from metainformant.ontology.visualization import visualization  # noqa: F401

    def test_ontology_core_go_functions(self) -> None:
        """Ontology core.go should expose load_go_obo."""
        from metainformant.ontology.core.go import load_go_obo  # noqa: F401

    def test_ontology_core_obo_functions(self) -> None:
        """Ontology core.obo should expose parse_obo."""
        from metainformant.ontology.core.obo import parse_obo  # noqa: F401

    def test_ontology_types(self) -> None:
        """Ontology types should be accessible."""
        from metainformant.ontology.core.types import Ontology, Term  # noqa: F401

    def test_visualization_canonical_paths(self) -> None:
        """Visualization canonical subpackage paths should work."""
        from metainformant.visualization.config import palettes, themes  # noqa: F401
        from metainformant.visualization.dashboards import composite, interactive  # noqa: F401

    def test_visualization_config_themes(self) -> None:
        """Visualization config.themes should expose theme functions."""
        from metainformant.visualization.config.themes import get_theme, apply_theme  # noqa: F401

    def test_visualization_config_palettes(self) -> None:
        """Visualization config.palettes should expose palette functions."""
        from metainformant.visualization.config.palettes import categorical, chromosome_palette  # noqa: F401

    def test_visualization_dashboards_composite(self) -> None:
        """Visualization dashboards.composite should expose dashboard functions."""
        from metainformant.visualization.dashboards.composite import multi_panel, qc_summary  # noqa: F401

    def test_visualization_dashboards_interactive(self) -> None:
        """Visualization dashboards.interactive should expose interactive plot functions."""
        from metainformant.visualization.dashboards.interactive import interactive_scatter  # noqa: F401

    def test_gwas_visualization_canonical_paths(self) -> None:
        """GWAS visualization canonical subpackage paths should work."""
        from metainformant.gwas.visualization import (  # noqa: F401
            genomic,
            interactive,
            population,
            statistical,
        )

    def test_gwas_visualization_population(self) -> None:
        """GWAS population visualization should work from canonical path."""
        from metainformant.gwas.visualization.population import admixture_plot, pca_scree_plot  # noqa: F401

    def test_gwas_visualization_genomic(self) -> None:
        """GWAS genomic visualization should work from canonical path."""
        from metainformant.gwas.visualization.genomic import genome_wide_ld_heatmap  # noqa: F401

    def test_gwas_visualization_statistical(self) -> None:
        """GWAS statistical visualization should work from canonical path."""
        from metainformant.gwas.visualization.statistical import power_plot  # noqa: F401

    def test_gwas_visualization_interactive(self) -> None:
        """GWAS interactive visualization should work from canonical path."""
        from metainformant.gwas.visualization.interactive import gwas_summary_panel  # noqa: F401

    def test_gwas_visualization_general(self) -> None:
        """GWAS general visualization should still work."""
        from metainformant.gwas.visualization import manhattan_plot, qq_plot  # noqa: F401

    def test_information_metrics_canonical_paths(self) -> None:
        """Information metrics canonical subpackage paths should work."""
        from metainformant.information.metrics.advanced import channel  # noqa: F401
        from metainformant.information.metrics.core import syntactic  # noqa: F401

    def test_information_metrics_core_syntactic(self) -> None:
        """Information metrics core.syntactic should expose entropy functions."""
        from metainformant.information.metrics.core.syntactic import (  # noqa: F401
            shannon_entropy,
            mutual_information,
        )

    def test_information_metrics_core_continuous(self) -> None:
        """Information metrics core.continuous should expose continuous functions."""
        from metainformant.information.metrics.core.continuous import differential_entropy  # noqa: F401

    def test_information_metrics_core_estimation(self) -> None:
        """Information metrics core.estimation should expose estimator functions."""
        from metainformant.information.metrics.core.estimation import entropy_estimator  # noqa: F401

    def test_information_metrics_advanced_channel(self) -> None:
        """Information metrics advanced.channel should expose channel functions."""
        from metainformant.information.metrics.advanced.channel import channel_capacity  # noqa: F401

    def test_information_metrics_advanced_geometry(self) -> None:
        """Information metrics advanced.geometry should expose geometry functions."""
        from metainformant.information.metrics.advanced.geometry import fisher_rao_distance  # noqa: F401

    def test_information_metrics_advanced_semantic(self) -> None:
        """Information metrics advanced.semantic should expose semantic functions."""
        from metainformant.information.metrics.advanced.semantic import semantic_similarity  # noqa: F401

    def test_information_metrics_advanced_hypothesis(self) -> None:
        """Information metrics advanced.hypothesis should expose test functions."""
        from metainformant.information.metrics.advanced.hypothesis import mi_permutation_test  # noqa: F401

    def test_information_metrics_advanced_decomposition(self) -> None:
        """Information metrics advanced.decomposition should expose PID functions."""
        from metainformant.information.metrics.advanced.decomposition import (  # noqa: F401
            partial_information_decomposition,
        )

    def test_math_population_genetics_coalescent(self) -> None:
        """Math population_genetics.coalescent should be the canonical path."""
        from metainformant.math.population_genetics.coalescent import (  # noqa: F401
            watterson_theta,
            tajimas_D,
        )

    def test_rna_engine_progress_tracker(self) -> None:
        """RNA engine.progress_tracker should be the canonical path."""
        from metainformant.rna.engine.progress_tracker import ProgressTracker  # noqa: F401
