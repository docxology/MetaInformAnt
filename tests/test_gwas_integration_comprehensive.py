"""Comprehensive integration tests for the full GWAS pipeline.

Exercises every module in the GWAS package end-to-end with synthetic data:
metadata, heritability, association, structure, correction, fine-mapping,
and all visualization sub-modules (phenotype, LD, composite, geography,
interactive, fine-mapping, config).

NO MOCKING -- all real implementations, real numpy/matplotlib computations.
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pytest

from metainformant.gwas.analysis.association import association_test_linear
from metainformant.gwas.analysis.correction import bonferroni_correction, fdr_correction
from metainformant.gwas.analysis.heritability import (
    estimate_heritability,
    heritability_bar_chart,
    partition_heritability_by_chromosome,
)
from metainformant.gwas.analysis.structure import compute_kinship_matrix, compute_pca
from metainformant.gwas.data.metadata import (
    get_geographic_coordinates,
    get_population_labels,
    load_sample_metadata,
    merge_metadata_with_phenotypes,
    validate_metadata,
)
from metainformant.gwas.visualization.config import THEMES, PlotStyle, apply_style, get_style
from metainformant.gwas.visualization.genomic.ld import compute_ld_decay, ld_decay_plot, ld_heatmap_region
from metainformant.gwas.visualization.interactive.composite import (
    gwas_summary_panel,
    population_structure_panel,
)
from metainformant.gwas.visualization.interactive.finemapping import compute_credible_set, credible_set_plot
from metainformant.gwas.visualization.interactive.interactive import interactive_manhattan
from metainformant.gwas.visualization.interactive.phenotype import (
    genotype_phenotype_boxplot,
    phenotype_correlation_matrix,
    phenotype_distribution,
    phenotype_pca_correlation,
    top_hits_genotype_phenotype,
)
from metainformant.gwas.visualization.population.geography import population_count_map, sample_map

# ---------------------------------------------------------------------------
# Helpers to generate synthetic data independent of the VCF fixture system
# ---------------------------------------------------------------------------

_POPULATIONS = ["CarnicaEU", "LigusticaIT", "MelliferaUK"]


def _make_sample_ids(n: int) -> List[str]:
    """Generate sample IDs: BEE_001 .. BEE_n."""
    return [f"BEE_{i:03d}" for i in range(1, n + 1)]


def _assign_populations(sample_ids: List[str]) -> Dict[str, str]:
    """Round-robin assign samples to 3 populations."""
    labels: Dict[str, str] = {}
    for i, sid in enumerate(sample_ids):
        labels[sid] = _POPULATIONS[i % len(_POPULATIONS)]
    return labels


def _make_metadata_tsv(
    path: Path,
    sample_ids: List[str],
    pop_labels: Dict[str, str],
) -> Path:
    """Write a metadata TSV with sample_id, population, latitude, longitude."""
    path.parent.mkdir(parents=True, exist_ok=True)

    # Assign realistic-ish coords per population
    pop_coords = {
        "CarnicaEU": (46.0, 14.5),
        "LigusticaIT": (43.0, 12.0),
        "MelliferaUK": (51.5, -1.2),
    }

    with open(path, "w") as fh:
        fh.write("sample_id\tpopulation\tlatitude\tlongitude\n")
        rng = random.Random(42)
        for sid in sample_ids:
            pop = pop_labels[sid]
            base_lat, base_lon = pop_coords[pop]
            lat = base_lat + rng.uniform(-0.5, 0.5)
            lon = base_lon + rng.uniform(-0.5, 0.5)
            fh.write(f"{sid}\t{pop}\t{lat:.4f}\t{lon:.4f}\n")

    return path


def _make_genotype_matrix(
    n_samples: int,
    n_variants: int,
    seed: int = 42,
) -> List[List[int]]:
    """Generate a genotype matrix (samples x variants) with values 0/1/2.

    Returns the matrix in *samples-major* layout: outer list = samples,
    inner list = genotypes across variants.
    """
    rng = np.random.default_rng(seed)
    # Allele frequencies per variant
    afs = rng.uniform(0.1, 0.5, size=n_variants)
    matrix: List[List[int]] = []
    for _ in range(n_samples):
        row: List[int] = []
        for af in afs:
            r = rng.random()
            if r < (1 - af) ** 2:
                row.append(0)
            elif r < (1 - af) ** 2 + 2 * af * (1 - af):
                row.append(1)
            else:
                row.append(2)
        matrix.append(row)
    return matrix


def _make_phenotypes(
    genotype_matrix: List[List[int]],
    causal_indices: List[int],
    effect_size: float = 2.0,
    seed: int = 42,
) -> List[float]:
    """Generate continuous phenotypes with known causal effects."""
    rng = np.random.default_rng(seed)
    n_samples = len(genotype_matrix)
    phenotypes: List[float] = []
    for s in range(n_samples):
        val = 50.0 + rng.normal(0, 3.0)
        for ci in causal_indices:
            val += genotype_matrix[s][ci] * effect_size
        phenotypes.append(float(val))
    return phenotypes


def _transpose_genotype_matrix(
    samples_x_variants: List[List[int]],
) -> List[List[int]]:
    """Transpose from (samples x variants) to (variants x samples)."""
    n_samples = len(samples_x_variants)
    n_variants = len(samples_x_variants[0]) if n_samples > 0 else 0
    return [[samples_x_variants[s][v] for s in range(n_samples)] for v in range(n_variants)]


def _make_positions(n_variants: int, chrom_size: int = 500_000) -> List[int]:
    """Generate evenly spaced variant positions."""
    spacing = chrom_size // max(n_variants, 1)
    return [1000 + i * spacing for i in range(n_variants)]


def _make_chromosomes(n_variants: int, n_chroms: int = 3) -> List[int]:
    """Assign variants to chromosomes round-robin."""
    return [(i % n_chroms) + 1 for i in range(n_variants)]


def _make_assoc_results(
    genotypes_by_variant: List[List[int]],
    phenotypes: List[float],
    positions: List[int],
    chromosomes: List[int],
) -> List[Dict[str, Any]]:
    """Run real linear association tests per variant."""
    from metainformant.gwas.analysis.association import association_test_linear

    results: List[Dict[str, Any]] = []
    for v_idx in range(len(genotypes_by_variant)):
        genos = genotypes_by_variant[v_idx]
        r = association_test_linear(genos, phenotypes)
        r["variant_id"] = f"rs{v_idx + 1}"
        r["chrom"] = chromosomes[v_idx]
        r["pos"] = positions[v_idx]
        r["position"] = positions[v_idx]
        r["chromosome"] = chromosomes[v_idx]
        results.append(r)
    return results


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

N_SAMPLES = 100
N_VARIANTS = 50
N_CHROMS = 3
CAUSAL_INDICES = [5, 15, 35]
EFFECT_SIZE = 3.0


@pytest.fixture
def sample_ids() -> List[str]:
    return _make_sample_ids(N_SAMPLES)


@pytest.fixture
def pop_labels(sample_ids: List[str]) -> Dict[str, str]:
    return _assign_populations(sample_ids)


@pytest.fixture
def genotype_matrix() -> List[List[int]]:
    """Samples x variants."""
    return _make_genotype_matrix(N_SAMPLES, N_VARIANTS, seed=42)


@pytest.fixture
def genotypes_by_variant(genotype_matrix: List[List[int]]) -> List[List[int]]:
    """Variants x samples."""
    return _transpose_genotype_matrix(genotype_matrix)


@pytest.fixture
def positions() -> List[int]:
    return _make_positions(N_VARIANTS)


@pytest.fixture
def chromosomes() -> List[int]:
    return _make_chromosomes(N_VARIANTS, N_CHROMS)


@pytest.fixture
def phenotypes(genotype_matrix: List[List[int]]) -> List[float]:
    return _make_phenotypes(genotype_matrix, CAUSAL_INDICES, EFFECT_SIZE, seed=42)


@pytest.fixture
def phenotypes_second_trait(genotype_matrix: List[List[int]]) -> List[float]:
    """A second independent trait for multi-trait tests."""
    return _make_phenotypes(genotype_matrix, [10, 25], 1.5, seed=99)


@pytest.fixture
def assoc_results(
    genotypes_by_variant: List[List[int]],
    phenotypes: List[float],
    positions: List[int],
    chromosomes: List[int],
) -> List[Dict[str, Any]]:
    return _make_assoc_results(genotypes_by_variant, phenotypes, positions, chromosomes)


@pytest.fixture
def metadata_dict(sample_ids: List[str], pop_labels: Dict[str, str]) -> Dict[str, Dict[str, Any]]:
    """In-memory metadata dict matching sample_map / geography expectations."""
    rng = random.Random(42)
    pop_coords = {
        "CarnicaEU": (46.0, 14.5),
        "LigusticaIT": (43.0, 12.0),
        "MelliferaUK": (51.5, -1.2),
    }
    meta: Dict[str, Dict[str, Any]] = {}
    for sid in sample_ids:
        pop = pop_labels[sid]
        base_lat, base_lon = pop_coords[pop]
        meta[sid] = {
            "population": pop,
            "latitude": base_lat + rng.uniform(-0.5, 0.5),
            "longitude": base_lon + rng.uniform(-0.5, 0.5),
        }
    return meta


# ---------------------------------------------------------------------------
# Test class: Module imports
# ---------------------------------------------------------------------------


class TestModuleImports:
    """Verify all GWAS sub-modules and key exports are importable."""

    def test_metadata_module(self) -> None:
        """Test metadata module imports."""
        assert callable(load_sample_metadata)
        assert callable(validate_metadata)
        assert callable(get_population_labels)
        assert callable(get_geographic_coordinates)
        assert callable(merge_metadata_with_phenotypes)

    def test_heritability_module(self) -> None:
        """Test heritability module imports."""
        assert callable(estimate_heritability)
        assert callable(partition_heritability_by_chromosome)
        assert callable(heritability_bar_chart)

    def test_visualization_phenotype_module(self) -> None:
        """Test visualization phenotype module imports."""
        assert callable(phenotype_distribution)
        assert callable(phenotype_correlation_matrix)
        assert callable(genotype_phenotype_boxplot)
        assert callable(top_hits_genotype_phenotype)
        assert callable(phenotype_pca_correlation)

    def test_visualization_ld_module(self) -> None:
        from metainformant.gwas.visualization.genomic.ld import (
            compute_ld_decay,
            ld_decay_plot,
            ld_heatmap_region,
        )

        assert callable(compute_ld_decay)
        assert callable(ld_decay_plot)
        assert callable(ld_heatmap_region)

    def test_visualization_composite_module(self) -> None:
        from metainformant.gwas.visualization.interactive.composite import (
            gwas_summary_panel,
            population_structure_panel,
            top_hit_detail_panel,
        )

        assert callable(gwas_summary_panel)
        assert callable(population_structure_panel)
        assert callable(top_hit_detail_panel)

    def test_visualization_geography_module(self) -> None:
        from metainformant.gwas.visualization.population.geography import (
            allele_frequency_map,
            population_count_map,
            sample_map,
        )

        assert callable(sample_map)
        assert callable(allele_frequency_map)
        assert callable(population_count_map)

    def test_visualization_interactive_module(self) -> None:
        from metainformant.gwas.visualization.interactive.interactive import (
            interactive_manhattan,
            interactive_pca,
            interactive_volcano,
        )

        assert callable(interactive_manhattan)
        assert callable(interactive_pca)
        assert callable(interactive_volcano)

    def test_visualization_finemapping_module(self) -> None:
        from metainformant.gwas.visualization.interactive.finemapping import (
            compute_credible_set,
            conditional_analysis_plot,
            credible_set_plot,
            pip_vs_ld_plot,
        )

        assert callable(compute_credible_set)
        assert callable(credible_set_plot)
        assert callable(conditional_analysis_plot)
        assert callable(pip_vs_ld_plot)

    def test_config_module(self) -> None:
        from metainformant.gwas.visualization.config import (
            THEMES,
            PlotStyle,
            apply_style,
            get_style,
            style_from_config,
        )

        assert callable(get_style)
        assert callable(apply_style)
        assert callable(style_from_config)
        assert isinstance(THEMES, dict)
        assert isinstance(PlotStyle(), PlotStyle)

    def test_top_level_gwas_exports(self) -> None:
        """Verify commonly used functions are importable from the top-level gwas package."""
        from metainformant.gwas import (
            THEMES,
            PlotStyle,
            apply_style,
            association_test_linear,
            bonferroni_correction,
            compute_credible_set,
            compute_kinship_matrix,
            compute_ld_decay,
            compute_pca,
            credible_set_plot,
            estimate_heritability,
            fdr_correction,
            genotype_phenotype_boxplot,
            get_population_labels,
            get_style,
            gwas_summary_panel,
            interactive_manhattan,
            ld_decay_plot,
            load_sample_metadata,
            phenotype_correlation_matrix,
            phenotype_distribution,
            population_count_map,
            population_structure_panel,
            sample_map,
            validate_metadata,
        )

        assert callable(association_test_linear)
        assert callable(estimate_heritability)
        assert callable(interactive_manhattan)


# ---------------------------------------------------------------------------
# Test class: PlotStyle / config themes
# ---------------------------------------------------------------------------


class TestConfigThemes:
    """Verify all PlotStyle themes are valid and apply_style works."""

    def test_publication_theme(self) -> None:
        from metainformant.gwas.visualization.config import THEMES, PlotStyle, get_style

        style = get_style("publication")
        assert isinstance(style, PlotStyle)
        assert style.dpi == 300
        assert style.theme == "publication"

    def test_presentation_theme(self) -> None:
        from metainformant.gwas.visualization.config import PlotStyle, get_style

        style = get_style("presentation")
        assert isinstance(style, PlotStyle)
        assert style.dpi == 150
        assert style.font_size == 14
        assert style.title_size == 18
        assert style.figsize == (14, 8)

    def test_poster_theme(self) -> None:
        from metainformant.gwas.visualization.config import PlotStyle, get_style

        style = get_style("poster")
        assert isinstance(style, PlotStyle)
        assert style.dpi == 300
        assert style.point_size == 12.0
        assert style.font_size == 18
        assert style.title_size == 24
        assert style.figsize == (20, 12)

    def test_all_themes_in_registry(self) -> None:
        from metainformant.gwas.visualization.config import THEMES

        assert "publication" in THEMES
        assert "presentation" in THEMES
        assert "poster" in THEMES
        assert len(THEMES) == 3

    def test_unknown_theme_fallback(self) -> None:
        from metainformant.gwas.visualization.config import get_style

        style = get_style("nonexistent_theme")
        assert style.theme == "publication"

    def test_get_style_overrides(self) -> None:
        from metainformant.gwas.visualization.config import get_style

        style = get_style("publication", dpi=72, alpha=0.5)
        assert style.dpi == 72
        assert style.alpha == 0.5

    def test_apply_style_does_not_raise(self) -> None:
        from metainformant.gwas.visualization.config import apply_style, get_style

        try:
            import matplotlib  # noqa: F401

            style = get_style("presentation")
            apply_style(style)
        except ImportError:
            pytest.skip("matplotlib not available")

    def test_style_from_config(self) -> None:
        from metainformant.gwas.visualization.config import PlotStyle, style_from_config

        cfg: Dict[str, Any] = {
            "output": {
                "visualization": {
                    "dpi": 72,
                    "font_size": 8,
                    "theme": "custom",
                }
            }
        }
        style = style_from_config(cfg)
        assert isinstance(style, PlotStyle)
        assert style.dpi == 72
        assert style.font_size == 8
        assert style.theme == "custom"


# ---------------------------------------------------------------------------
# Test class: Full pipeline with synthetic data
# ---------------------------------------------------------------------------


class TestFullPipelineSynthetic:
    """End-to-end pipeline test exercising every new module with synthetic data."""

    # -- Metadata ----------------------------------------------------------

    def test_metadata_load_and_validate(
        self,
        tmp_path: Path,
        sample_ids: List[str],
        pop_labels: Dict[str, str],
    ) -> None:
        """Load metadata TSV, validate against sample list, extract labels & coords."""
        from metainformant.gwas.data.metadata import (
            get_geographic_coordinates,
            get_population_labels,
            load_sample_metadata,
            validate_metadata,
        )

        tsv_path = _make_metadata_tsv(tmp_path / "meta.tsv", sample_ids, pop_labels)

        # Load
        meta_result = load_sample_metadata(tsv_path)
        assert meta_result["status"] == "success"
        assert meta_result["n_samples"] == N_SAMPLES
        assert "population" in meta_result["columns"]

        # Validate
        val = validate_metadata(meta_result["metadata"], sample_ids)
        assert val["status"] == "success"
        assert val["valid"] is True
        assert val["completeness"] == 1.0

        # Population labels
        labels = get_population_labels(meta_result["metadata"])
        assert len(labels) == N_SAMPLES
        assert set(labels.values()) == set(_POPULATIONS)

        # Coordinates
        coords = get_geographic_coordinates(meta_result["metadata"])
        assert len(coords) == N_SAMPLES
        assert all("latitude" in c and "longitude" in c for c in coords)

    # -- Association testing -----------------------------------------------

    def test_association_testing(
        self,
        genotypes_by_variant: List[List[int]],
        phenotypes: List[float],
    ) -> None:
        """Run linear association on every variant and verify real results."""
        from metainformant.gwas.analysis.association import association_test_linear

        p_values: List[float] = []
        for genos in genotypes_by_variant:
            r = association_test_linear(genos, phenotypes)
            assert r["status"] == "success"
            assert 0.0 <= r["p_value"] <= 1.0
            p_values.append(r["p_value"])

        # p-values should not all be identical (real computation)
        assert len(set(round(p, 8) for p in p_values)) > 1

    # -- Multiple testing correction ---------------------------------------

    def test_multiple_testing_correction(self, assoc_results: List[Dict[str, Any]]) -> None:
        from metainformant.gwas.analysis.correction import bonferroni_correction, fdr_correction

        p_values = [r["p_value"] for r in assoc_results if r.get("status") == "success"]
        assert len(p_values) == N_VARIANTS

        bonf = bonferroni_correction(p_values)
        assert bonf["status"] == "success"
        assert bonf["n_tests"] == N_VARIANTS
        assert bonf["corrected_alpha"] == pytest.approx(0.05 / N_VARIANTS)

        fdr = fdr_correction(p_values)
        assert fdr["status"] == "success"
        assert fdr["n_tests"] == N_VARIANTS

    # -- Population structure (PCA + kinship) ------------------------------

    def test_pca_and_kinship(self, genotype_matrix: List[List[int]]) -> None:
        from metainformant.gwas.analysis.structure import compute_kinship_matrix, compute_pca

        pca_result = compute_pca(genotype_matrix, n_components=5)
        assert pca_result["status"] == "success"
        assert len(pca_result["pcs"]) == N_SAMPLES
        assert len(pca_result["pcs"][0]) == 5
        assert len(pca_result["explained_variance_ratio"]) == 5

        kinship_result = compute_kinship_matrix(genotype_matrix)
        assert kinship_result["status"] == "success"
        km = kinship_result["kinship_matrix"]
        assert len(km) == N_SAMPLES
        assert len(km[0]) == N_SAMPLES

    # -- Heritability estimation -------------------------------------------

    def test_heritability_estimation(
        self,
        genotype_matrix: List[List[int]],
        phenotypes: List[float],
    ) -> None:
        from metainformant.gwas.analysis.heritability import (
            estimate_heritability,
            partition_heritability_by_chromosome,
        )
        from metainformant.gwas.analysis.structure import compute_kinship_matrix

        kinship_result = compute_kinship_matrix(genotype_matrix)
        km = np.array(kinship_result["kinship_matrix"])

        h2_result = estimate_heritability(km, phenotypes)
        assert h2_result["status"] == "success"
        assert 0.0 <= h2_result["h2"] <= 1.0
        assert h2_result["h2_se"] >= 0.0
        assert h2_result["n_samples"] == N_SAMPLES

    def test_heritability_partitioning(
        self,
        genotype_matrix: List[List[int]],
        phenotypes: List[float],
        chromosomes: List[int],
    ) -> None:
        from metainformant.gwas.analysis.heritability import partition_heritability_by_chromosome
        from metainformant.gwas.analysis.structure import compute_kinship_matrix

        # Build per-chromosome genotype sub-matrices and kinship matrices
        unique_chroms = sorted(set(chromosomes))
        chrom_kinship: Dict[int, Any] = {}
        for chrom in unique_chroms:
            variant_indices = [v for v in range(N_VARIANTS) if chromosomes[v] == chrom]
            sub_matrix = [[genotype_matrix[s][v] for v in variant_indices] for s in range(N_SAMPLES)]
            kr = compute_kinship_matrix(sub_matrix)
            if kr["status"] == "success":
                chrom_kinship[chrom] = np.array(kr["kinship_matrix"])

        part_result = partition_heritability_by_chromosome(chrom_kinship, phenotypes)
        assert part_result["status"] == "success"
        assert part_result["n_chromosomes"] == len(unique_chroms)
        assert 0.0 <= part_result["total_h2"] <= 1.0

    # -- Fine-mapping (credible sets) --------------------------------------

    def test_credible_set_computation(self, assoc_results: List[Dict[str, Any]]) -> None:
        from metainformant.gwas.visualization.interactive.finemapping import compute_credible_set

        cs = compute_credible_set(assoc_results, credible_level=0.95)
        assert cs["status"] == "success"
        assert len(cs["pips"]) == N_VARIANTS
        assert cs["credible_set_size"] > 0
        assert cs["credible_set_size"] <= N_VARIANTS
        assert cs["cumulative_probability"] >= 0.95

        # PIPs should sum to ~1.0
        assert abs(sum(cs["pips"]) - 1.0) < 1e-6

    # -- Visualization: phenotype distribution -----------------------------

    def test_phenotype_distribution_plot(
        self,
        tmp_path: Path,
        phenotypes: List[float],
        sample_ids: List[str],
        pop_labels: Dict[str, str],
    ) -> None:
        from metainformant.gwas.visualization.interactive.phenotype import phenotype_distribution

        out = tmp_path / "pheno_dist.png"
        result = phenotype_distribution(
            phenotype_values=phenotypes,
            trait_name="Varroa Resistance",
            output_file=out,
            by_population=pop_labels,
            sample_ids=sample_ids,
        )
        assert result["status"] == "success"
        assert result["n_samples"] == N_SAMPLES
        assert out.exists()
        assert out.stat().st_size > 0

    # -- Visualization: phenotype correlation matrix -----------------------

    def test_phenotype_correlation_matrix(
        self,
        tmp_path: Path,
        phenotypes: List[float],
        phenotypes_second_trait: List[float],
    ) -> None:
        from metainformant.gwas.visualization.interactive.phenotype import phenotype_correlation_matrix

        out = tmp_path / "corr_matrix.png"
        traits = {
            "Varroa Resistance": phenotypes,
            "Honey Production": phenotypes_second_trait,
        }
        result = phenotype_correlation_matrix(traits, output_file=out)
        assert result["status"] == "success"
        assert result["n_traits"] == 2
        assert "correlation_matrix" in result
        assert out.exists()

    # -- Visualization: genotype-phenotype boxplot -------------------------

    def test_genotype_phenotype_boxplot(
        self,
        tmp_path: Path,
        genotypes_by_variant: List[List[int]],
        phenotypes: List[float],
        assoc_results: List[Dict[str, Any]],
    ) -> None:
        from metainformant.gwas.visualization.interactive.phenotype import genotype_phenotype_boxplot

        # Find the top hit by p-value
        sorted_results = sorted(
            enumerate(assoc_results),
            key=lambda x: x[1].get("p_value", 1.0),
        )
        top_idx, top_result = sorted_results[0]

        out = tmp_path / "geno_pheno_box.png"
        result = genotype_phenotype_boxplot(
            genotypes=genotypes_by_variant[top_idx],
            phenotypes=phenotypes,
            variant_id=top_result.get("variant_id", f"rs{top_idx}"),
            output_file=out,
        )
        assert result["status"] == "success"
        assert result["n_samples"] == N_SAMPLES
        assert out.exists()

    # -- Visualization: LD decay -------------------------------------------

    def test_ld_decay_plot(
        self,
        tmp_path: Path,
        genotypes_by_variant: List[List[int]],
        positions: List[int],
        chromosomes: List[int],
    ) -> None:
        from metainformant.gwas.visualization.genomic.ld import compute_ld_decay, ld_decay_plot

        decay = compute_ld_decay(
            genotypes_by_variant=genotypes_by_variant,
            positions=positions,
            max_distance=200_000,
            n_bins=20,
            chromosomes=chromosomes,
        )
        assert decay["status"] == "success"
        assert decay["total_pairs"] > 0
        assert len(decay["bin_centers"]) == 20

        out = tmp_path / "ld_decay.png"
        plot_result = ld_decay_plot(decay, output_file=out, title="LD Decay (Test)")
        assert plot_result["status"] == "success"
        assert out.exists()

    # -- Visualization: LD heatmap -----------------------------------------

    def test_ld_heatmap_region(
        self,
        tmp_path: Path,
        genotypes_by_variant: List[List[int]],
        positions: List[int],
    ) -> None:
        from metainformant.gwas.visualization.genomic.ld import ld_heatmap_region

        # Use first 15 variants for a manageable heatmap
        out = tmp_path / "ld_heatmap.png"
        result = ld_heatmap_region(
            genotypes_by_variant=genotypes_by_variant[:15],
            positions=positions[:15],
            output_file=out,
            title="LD Heatmap (Test Region)",
        )
        assert result["status"] == "success"
        assert result["n_variants"] == 15
        assert out.exists()

    # -- Visualization: PCA multi-panel ------------------------------------

    def test_pca_multi_panel(
        self,
        tmp_path: Path,
        genotype_matrix: List[List[int]],
        sample_ids: List[str],
        pop_labels: Dict[str, str],
    ) -> None:
        from metainformant.gwas.analysis.structure import compute_pca
        from metainformant.gwas.visualization.population.population_pca import pca_multi_panel

        pca_result = compute_pca(genotype_matrix, n_components=5)
        assert pca_result["status"] == "success"

        # Build metadata dict for coloring
        meta = {sid: {"population": pop_labels[sid]} for sid in sample_ids}

        out = tmp_path / "pca_multi.png"
        result = pca_multi_panel(
            pca_data={
                "pcs": pca_result["pcs"],
                "explained_variance_ratio": pca_result["explained_variance_ratio"],
                "sample_ids": sample_ids,
            },
            metadata=meta,
            output_file=out,
            title="PCA Multi-Panel (Test)",
        )
        assert result["status"] == "success"
        assert out.exists()

    # -- Visualization: GWAS summary panel (composite) ---------------------

    def test_gwas_summary_panel(
        self,
        tmp_path: Path,
        assoc_results: List[Dict[str, Any]],
        genotype_matrix: List[List[int]],
    ) -> None:
        from metainformant.gwas.analysis.structure import compute_kinship_matrix, compute_pca
        from metainformant.gwas.visualization.interactive.composite import gwas_summary_panel

        pca_result = compute_pca(genotype_matrix, n_components=3)
        kinship_result = compute_kinship_matrix(genotype_matrix)

        out = tmp_path / "gwas_summary.png"
        result = gwas_summary_panel(
            assoc_results=assoc_results,
            pca_data={
                "pcs": pca_result["pcs"],
                "explained_variance_ratio": pca_result["explained_variance_ratio"],
            },
            kinship_matrix=kinship_result["kinship_matrix"],
            output_file=out,
            title="GWAS Summary (Test)",
        )
        assert result["status"] == "success"
        assert result["panels_generated"] == 4
        assert out.exists()

    # -- Visualization: population structure panel -------------------------

    def test_population_structure_panel(
        self,
        tmp_path: Path,
        genotype_matrix: List[List[int]],
        sample_ids: List[str],
        pop_labels: Dict[str, str],
    ) -> None:
        from metainformant.gwas.analysis.structure import compute_kinship_matrix, compute_pca
        from metainformant.gwas.visualization.interactive.composite import population_structure_panel

        pca_result = compute_pca(genotype_matrix, n_components=5)
        kinship_result = compute_kinship_matrix(genotype_matrix)
        meta = {sid: {"population": pop_labels[sid]} for sid in sample_ids}

        out = tmp_path / "pop_structure.png"
        result = population_structure_panel(
            pca_data={
                "pcs": pca_result["pcs"],
                "explained_variance_ratio": pca_result["explained_variance_ratio"],
                "sample_ids": sample_ids,
            },
            kinship_matrix=kinship_result["kinship_matrix"],
            metadata=meta,
            output_file=out,
            title="Population Structure (Test)",
        )
        assert result["status"] == "success"
        assert result["panels_generated"] == 4
        assert out.exists()

    # -- Visualization: interactive Manhattan (HTML) -----------------------

    def test_interactive_manhattan(
        self,
        tmp_path: Path,
        assoc_results: List[Dict[str, Any]],
    ) -> None:
        from metainformant.gwas.visualization.interactive.interactive import interactive_manhattan

        out = tmp_path / "manhattan.html"
        result = interactive_manhattan(
            assoc_results=assoc_results,
            output_file=out,
            title="Interactive Manhattan (Test)",
        )
        assert result["status"] == "success"
        assert result["n_variants"] == N_VARIANTS
        assert out.exists()
        content = out.read_text()
        assert len(content) > 100

    # -- Visualization: interactive volcano --------------------------------

    def test_interactive_volcano(
        self,
        tmp_path: Path,
        assoc_results: List[Dict[str, Any]],
    ) -> None:
        from metainformant.gwas.visualization.interactive.interactive import interactive_volcano

        out = tmp_path / "volcano.html"
        result = interactive_volcano(
            assoc_results=assoc_results,
            output_file=out,
            title="Interactive Volcano (Test)",
        )
        assert result["status"] == "success"
        assert result["n_variants"] == N_VARIANTS
        assert out.exists()

    # -- Visualization: interactive PCA ------------------------------------

    def test_interactive_pca(
        self,
        tmp_path: Path,
        genotype_matrix: List[List[int]],
        sample_ids: List[str],
        pop_labels: Dict[str, str],
    ) -> None:
        from metainformant.gwas.analysis.structure import compute_pca
        from metainformant.gwas.visualization.interactive.interactive import interactive_pca

        pca_result = compute_pca(genotype_matrix, n_components=5)
        meta = {sid: {"population": pop_labels[sid]} for sid in sample_ids}

        out = tmp_path / "pca_3d.html"
        result = interactive_pca(
            pca_data={
                "pcs": pca_result["pcs"],
                "explained_variance_ratio": pca_result["explained_variance_ratio"],
                "sample_ids": sample_ids,
            },
            metadata=meta,
            output_file=out,
            title="Interactive 3D PCA (Test)",
        )
        assert result["status"] == "success"
        assert result["n_samples"] == N_SAMPLES
        assert out.exists()

    # -- Visualization: sample geography map -------------------------------

    def test_sample_map(
        self,
        tmp_path: Path,
        metadata_dict: Dict[str, Dict[str, Any]],
    ) -> None:
        from metainformant.gwas.visualization.population.geography import sample_map

        out = tmp_path / "sample_map.png"
        result = sample_map(
            metadata=metadata_dict,
            output_file=out,
            color_by="population",
            title="Sample Map (Test)",
        )
        assert result["status"] == "success"
        assert result["n_samples_plotted"] == N_SAMPLES
        assert result["n_skipped"] == 0
        assert out.exists()

    # -- Visualization: population count map -------------------------------

    def test_population_count_map(
        self,
        tmp_path: Path,
        metadata_dict: Dict[str, Dict[str, Any]],
    ) -> None:
        from metainformant.gwas.visualization.population.geography import population_count_map

        out = tmp_path / "pop_count_map.png"
        result = population_count_map(
            metadata=metadata_dict,
            output_file=out,
            title="Population Count Map (Test)",
        )
        assert result["status"] == "success"
        assert result["n_populations"] == len(_POPULATIONS)
        assert out.exists()

    # -- Visualization: credible set plot ----------------------------------

    def test_credible_set_plot(
        self,
        tmp_path: Path,
        assoc_results: List[Dict[str, Any]],
    ) -> None:
        from metainformant.gwas.visualization.interactive.finemapping import credible_set_plot

        out = tmp_path / "credible_set.png"
        result = credible_set_plot(
            assoc_results=assoc_results,
            output_file=out,
            credible_level=0.95,
            title="Credible Set (Test)",
        )
        assert result["status"] == "success"
        assert result["credible_set_size"] > 0
        assert out.exists()

    # -- Visualization: PIP vs LD ------------------------------------------

    def test_pip_vs_ld_plot(
        self,
        tmp_path: Path,
        assoc_results: List[Dict[str, Any]],
    ) -> None:
        from metainformant.gwas.visualization.interactive.finemapping import compute_credible_set, pip_vs_ld_plot

        cs = compute_credible_set(assoc_results)
        pips = cs["pips"]

        # Synthesize LD values with "lead" variant (highest PIP)
        rng = np.random.default_rng(42)
        ld_with_lead = rng.uniform(0.0, 1.0, size=N_VARIANTS).tolist()
        # Lead variant has LD=1.0 with itself
        lead_idx = int(np.argmax(pips))
        ld_with_lead[lead_idx] = 1.0

        out = tmp_path / "pip_vs_ld.png"
        result = pip_vs_ld_plot(
            pips=pips,
            ld_with_lead=ld_with_lead,
            output_file=out,
            variant_ids=[f"rs{i + 1}" for i in range(N_VARIANTS)],
            title="PIP vs LD (Test)",
        )
        assert result["status"] == "success"
        assert result["n_variants"] == N_VARIANTS
        assert out.exists()

    # -- Visualization: heritability bar chart ------------------------------

    def test_heritability_bar_chart(
        self,
        tmp_path: Path,
        genotype_matrix: List[List[int]],
        phenotypes: List[float],
        chromosomes: List[int],
    ) -> None:
        from metainformant.gwas.analysis.heritability import (
            heritability_bar_chart,
            partition_heritability_by_chromosome,
        )
        from metainformant.gwas.analysis.structure import compute_kinship_matrix

        unique_chroms = sorted(set(chromosomes))
        chrom_kinship: Dict[int, Any] = {}
        for chrom in unique_chroms:
            variant_indices = [v for v in range(N_VARIANTS) if chromosomes[v] == chrom]
            sub_matrix = [[genotype_matrix[s][v] for v in variant_indices] for s in range(N_SAMPLES)]
            kr = compute_kinship_matrix(sub_matrix)
            if kr["status"] == "success":
                chrom_kinship[chrom] = np.array(kr["kinship_matrix"])

        h2_data = partition_heritability_by_chromosome(chrom_kinship, phenotypes)
        assert h2_data["status"] == "success"

        out = tmp_path / "h2_bar.png"
        result = heritability_bar_chart(h2_data, output_file=out)
        assert result["status"] == "success"
        assert out.exists()

    # -- Visualization: phenotype PCA correlation --------------------------

    def test_phenotype_pca_correlation(
        self,
        tmp_path: Path,
        genotype_matrix: List[List[int]],
        phenotypes: List[float],
    ) -> None:
        from metainformant.gwas.analysis.structure import compute_pca
        from metainformant.gwas.visualization.interactive.phenotype import phenotype_pca_correlation

        pca_result = compute_pca(genotype_matrix, n_components=3)
        assert pca_result["status"] == "success"

        out = tmp_path / "pheno_pca_corr.png"
        result = phenotype_pca_correlation(
            pcs=pca_result["pcs"],
            phenotype_values=phenotypes,
            output_file=out,
            trait_name="Varroa Resistance",
        )
        assert result["status"] == "success"
        assert result["n_samples"] == N_SAMPLES
        assert -1.0 <= result["r_pc1"] <= 1.0
        assert -1.0 <= result["r_pc2"] <= 1.0
        assert out.exists()

    # -- PlotStyle application works in pipeline context -------------------

    def test_apply_plot_style(self) -> None:
        """Applying a style should not error and should set rcParams."""
        from metainformant.gwas.visualization.config import apply_style, get_style

        try:
            import matplotlib as mpl

            style = get_style("poster")
            apply_style(style)
            assert mpl.rcParams["font.size"] == 18
            assert mpl.rcParams["axes.titlesize"] == 24
        except ImportError:
            pytest.skip("matplotlib not available")

    # -- Full pipeline orchestration test ----------------------------------

    def test_full_pipeline_synthetic(
        self,
        tmp_path: Path,
        sample_ids: List[str],
        pop_labels: Dict[str, str],
        genotype_matrix: List[List[int]],
        genotypes_by_variant: List[List[int]],
        phenotypes: List[float],
        phenotypes_second_trait: List[float],
        positions: List[int],
        chromosomes: List[int],
        metadata_dict: Dict[str, Dict[str, Any]],
    ) -> None:
        """Run the full pipeline end-to-end and verify all outputs exist."""
        out_dir = tmp_path / "full_pipeline"
        out_dir.mkdir(parents=True, exist_ok=True)

        # 1. Create metadata TSV & load it
        from metainformant.gwas.data.metadata import load_sample_metadata, validate_metadata

        tsv_path = _make_metadata_tsv(out_dir / "metadata.tsv", sample_ids, pop_labels)
        meta_result = load_sample_metadata(tsv_path)
        assert meta_result["status"] == "success"

        val = validate_metadata(meta_result["metadata"], sample_ids)
        assert val["valid"] is True

        # 2. Association testing
        from metainformant.gwas.analysis.association import association_test_linear

        assoc_results: List[Dict[str, Any]] = []
        for v_idx in range(N_VARIANTS):
            r = association_test_linear(genotypes_by_variant[v_idx], phenotypes)
            r["variant_id"] = f"rs{v_idx + 1}"
            r["chrom"] = chromosomes[v_idx]
            r["pos"] = positions[v_idx]
            r["position"] = positions[v_idx]
            r["chromosome"] = chromosomes[v_idx]
            assoc_results.append(r)

        p_values = [r["p_value"] for r in assoc_results]
        assert all(0.0 <= p <= 1.0 for p in p_values)

        # 3. Multiple testing correction
        from metainformant.gwas.analysis.correction import bonferroni_correction

        bonf = bonferroni_correction(p_values)
        assert bonf["status"] == "success"

        # 4. PCA & kinship
        from metainformant.gwas.analysis.structure import compute_kinship_matrix, compute_pca

        pca_result = compute_pca(genotype_matrix, n_components=5)
        assert pca_result["status"] == "success"

        kinship_result = compute_kinship_matrix(genotype_matrix)
        assert kinship_result["status"] == "success"
        km = np.array(kinship_result["kinship_matrix"])

        # 5. Heritability
        from metainformant.gwas.analysis.heritability import estimate_heritability

        h2 = estimate_heritability(km, phenotypes)
        assert h2["status"] == "success"
        assert 0.0 <= h2["h2"] <= 1.0

        # 6. Credible sets (fine-mapping)
        from metainformant.gwas.visualization.interactive.finemapping import compute_credible_set

        cs = compute_credible_set(assoc_results)
        assert cs["status"] == "success"
        assert cs["credible_set_size"] > 0

        # 7. Phenotype distribution
        from metainformant.gwas.visualization.interactive.phenotype import phenotype_distribution

        pheno_dist_path = out_dir / "pheno_dist.png"
        pd_result = phenotype_distribution(phenotypes, output_file=pheno_dist_path, trait_name="Varroa")
        assert pd_result["status"] == "success"
        assert pheno_dist_path.exists()

        # 8. Correlation matrix
        from metainformant.gwas.visualization.interactive.phenotype import phenotype_correlation_matrix

        corr_path = out_dir / "corr.png"
        corr_result = phenotype_correlation_matrix(
            {"Varroa": phenotypes, "Honey": phenotypes_second_trait},
            output_file=corr_path,
        )
        assert corr_result["status"] == "success"
        assert corr_path.exists()

        # 9. Genotype-phenotype boxplot for top hit
        from metainformant.gwas.visualization.interactive.phenotype import genotype_phenotype_boxplot

        sorted_assoc = sorted(enumerate(assoc_results), key=lambda x: x[1].get("p_value", 1.0))
        top_idx = sorted_assoc[0][0]
        box_path = out_dir / "geno_pheno.png"
        box_result = genotype_phenotype_boxplot(
            genotypes_by_variant[top_idx],
            phenotypes,
            variant_id=f"rs{top_idx + 1}",
            output_file=box_path,
        )
        assert box_result["status"] == "success"
        assert box_path.exists()

        # 10. LD decay
        from metainformant.gwas.visualization.genomic.ld import compute_ld_decay, ld_decay_plot

        decay = compute_ld_decay(genotypes_by_variant, positions, max_distance=200_000, n_bins=20)
        ld_path = out_dir / "ld_decay.png"
        ld_result = ld_decay_plot(decay, output_file=ld_path)
        assert ld_result["status"] == "success"
        assert ld_path.exists()

        # 11. PCA multi-panel
        from metainformant.gwas.visualization.population.population_pca import pca_multi_panel

        pca_panel_path = out_dir / "pca_multi.png"
        meta_for_panel = {sid: {"population": pop_labels[sid]} for sid in sample_ids}
        pca_panel_result = pca_multi_panel(
            pca_data={
                "pcs": pca_result["pcs"],
                "explained_variance_ratio": pca_result["explained_variance_ratio"],
                "sample_ids": sample_ids,
            },
            metadata=meta_for_panel,
            output_file=pca_panel_path,
        )
        assert pca_panel_result["status"] == "success"
        assert pca_panel_path.exists()

        # 12. GWAS summary panel (composite)
        from metainformant.gwas.visualization.interactive.composite import gwas_summary_panel

        summary_path = out_dir / "summary.png"
        summary_result = gwas_summary_panel(
            assoc_results=assoc_results,
            pca_data={
                "pcs": pca_result["pcs"],
                "explained_variance_ratio": pca_result["explained_variance_ratio"],
            },
            kinship_matrix=km,
            output_file=summary_path,
        )
        assert summary_result["status"] == "success"
        assert summary_path.exists()

        # 13. Interactive Manhattan (HTML)
        from metainformant.gwas.visualization.interactive.interactive import interactive_manhattan

        manhattan_path = out_dir / "manhattan.html"
        manhattan_result = interactive_manhattan(assoc_results, output_file=manhattan_path)
        assert manhattan_result["status"] == "success"
        assert manhattan_path.exists()

        # 14. Sample geography map
        from metainformant.gwas.visualization.population.geography import sample_map

        map_path = out_dir / "sample_map.png"
        map_result = sample_map(metadata=metadata_dict, output_file=map_path)
        assert map_result["status"] == "success"
        assert map_path.exists()

        # 15. Apply PlotStyle
        from metainformant.gwas.visualization.config import apply_style, get_style

        try:
            import matplotlib  # noqa: F401

            style = get_style("publication")
            apply_style(style)
        except ImportError:
            pass  # matplotlib optional for apply_style test

        # Verify all output files exist
        expected_files = [
            pheno_dist_path,
            corr_path,
            box_path,
            ld_path,
            pca_panel_path,
            summary_path,
            manhattan_path,
            map_path,
        ]
        for fp in expected_files:
            assert fp.exists(), f"Expected output file missing: {fp}"
            assert fp.stat().st_size > 0, f"Output file is empty: {fp}"
