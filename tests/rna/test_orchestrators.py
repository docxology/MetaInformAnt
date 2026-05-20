"""Integration tests for WithinSpeciesOrchestrator and AcrossSpeciesOrchestrator.

Validates end-to-end functionality of both orchestrators using synthetic data,
ensuring all analysis methods (PCA, DE, ortholog mapping, conservation scoring,
divergence matrix computation) work correctly from orchestrator entry-points.
"""

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis.across_species_orchestrator import AcrossSpeciesOrchestrator
from metainformant.rna.analysis.within_species_orchestrator import WithinSpeciesOrchestrator

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def synthetic_abundance(tmp_path):
    """Create a synthetic abundance matrix TSV file."""
    np.random.seed(42)
    genes = [f"gene_{i}" for i in range(100)]
    samples = [f"SRR{1000 + i}" for i in range(12)]
    data = np.random.poisson(lam=50, size=(100, 12)).astype(float)
    # Add clear differential signal for first 10 genes
    data[:10, 6:] *= 3  # Higher in second group
    df = pd.DataFrame(data, index=genes, columns=samples)
    path = tmp_path / "abundance.tsv"
    df.to_csv(path, sep="\t")
    return path, df


@pytest.fixture
def synthetic_metadata(tmp_path):
    """Create a synthetic metadata TSV file with tissue and caste columns."""
    samples = [f"SRR{1000 + i}" for i in range(12)]
    meta = pd.DataFrame(
        {
            "run": samples,
            "tissue": ["brain"] * 6 + ["antenna"] * 6,
            "caste": ["worker"] * 4 + ["queen"] * 4 + ["worker"] * 2 + ["queen"] * 2,
            "sex": ["female"] * 12,
        }
    )
    path = tmp_path / "metadata.tsv"
    meta.to_csv(path, sep="\t", index=False)
    return path, meta


@pytest.fixture
def three_species_data(tmp_path):
    """Create synthetic expression data for 3 species with shared orthologs."""
    np.random.seed(123)
    shared_genes = [f"ortho_{i}" for i in range(50)]

    species_data = {}
    for sp_name in ["Apis_mellifera", "Pogonomyrmex_barbatus", "Solenopsis_invicta"]:
        samples = [f"{sp_name}_s{i}" for i in range(5)]
        data = np.random.lognormal(mean=3, sigma=1.5, size=(50, 5))
        df = pd.DataFrame(data, index=shared_genes, columns=samples)
        path = tmp_path / f"{sp_name}_abundance.tsv"
        df.to_csv(path, sep="\t")
        species_data[sp_name] = (path, df)

    return species_data


@pytest.fixture
def ortholog_table(tmp_path, three_species_data):
    """Create a synthetic ortholog table mapping genes across 3 species."""
    shared_genes = [f"ortho_{i}" for i in range(50)]

    # Create an OrthoFinder-style table
    orth_df = pd.DataFrame(
        {
            "Orthogroup": [f"OG{i:04d}" for i in range(50)],
            "Apis_mellifera": shared_genes,
            "Pogonomyrmex_barbatus": shared_genes,
            "Solenopsis_invicta": shared_genes,
        }
    )
    path = tmp_path / "orthologs.tsv"
    orth_df.to_csv(path, sep="\t", index=False)
    return path, orth_df


# =============================================================================
# WithinSpeciesOrchestrator Tests
# =============================================================================


class TestWithinSpeciesOrchestrator:
    """Test suite for WithinSpeciesOrchestrator."""

    def test_instantiation(self, tmp_path, synthetic_abundance, synthetic_metadata):
        """Test that orchestrator can be instantiated with valid paths."""
        abundance_path, _ = synthetic_abundance
        metadata_path, _ = synthetic_metadata
        output_dir = tmp_path / "output"

        orch = WithinSpeciesOrchestrator(
            species_name="test_species",
            abundance_path=abundance_path,
            metadata_path=metadata_path,
            output_dir=output_dir,
        )
        assert orch.species_name == "test_species"
        assert output_dir.exists()

    def test_load_data(self, tmp_path, synthetic_abundance, synthetic_metadata):
        """Test data loading aligns samples correctly."""
        abundance_path, _ = synthetic_abundance
        metadata_path, _ = synthetic_metadata
        output_dir = tmp_path / "output"

        orch = WithinSpeciesOrchestrator("test_sp", abundance_path, metadata_path, output_dir)
        orch.load_data()

        assert orch.counts_df is not None
        assert orch.metadata_df is not None
        assert orch.counts_df.shape[1] == len(orch.metadata_df)
        # All column names in counts should be in metadata index
        assert set(orch.counts_df.columns) == set(orch.metadata_df.index)

    def test_run_pca(self, tmp_path, synthetic_abundance, synthetic_metadata):
        """Test PCA produces valid output with correct dimensionality."""
        abundance_path, _ = synthetic_abundance
        metadata_path, _ = synthetic_metadata
        output_dir = tmp_path / "output"

        orch = WithinSpeciesOrchestrator("test_sp", abundance_path, metadata_path, output_dir)
        orch.load_data()
        pca_res = orch.run_pca(n_components=2)

        assert "transformed" in pca_res
        assert "loadings" in pca_res
        assert "explained_variance_ratio" in pca_res
        assert pca_res["transformed"].shape[1] == 2
        assert pca_res["transformed"].shape[0] == 12  # number of samples

        # Check output file was written
        out_file = output_dir / "test_sp_pca_coordinates.tsv"
        assert out_file.exists()

    def test_run_differential_expression_tissue(self, tmp_path, synthetic_abundance, synthetic_metadata):
        """Test DE analysis with tissue condition (2 groups: brain vs antenna)."""
        abundance_path, _ = synthetic_abundance
        metadata_path, _ = synthetic_metadata
        output_dir = tmp_path / "output"

        orch = WithinSpeciesOrchestrator("test_sp", abundance_path, metadata_path, output_dir)
        orch.load_data()
        de_result = orch.run_differential_expression("tissue")

        assert de_result is not None
        assert "log2_fold_change" in de_result.columns
        assert "regulation" in de_result.columns  # from prepare_volcano_data
        assert len(de_result) > 0

        # Check output file was written
        out_file = output_dir / "test_sp_DE_tissue.tsv"
        assert out_file.exists()

    def test_de_skips_single_condition(self, tmp_path, synthetic_abundance, synthetic_metadata):
        """Test DE gracefully returns None for a column with only 1 unique value."""
        abundance_path, _ = synthetic_abundance
        metadata_path, _ = synthetic_metadata
        output_dir = tmp_path / "output"

        orch = WithinSpeciesOrchestrator("test_sp", abundance_path, metadata_path, output_dir)
        orch.load_data()
        # 'sex' column has only "female" — should return None
        result = orch.run_differential_expression("sex")
        assert result is None

    def test_de_skips_missing_column(self, tmp_path, synthetic_abundance, synthetic_metadata):
        """Test DE gracefully returns None for non-existent condition column."""
        abundance_path, _ = synthetic_abundance
        metadata_path, _ = synthetic_metadata
        output_dir = tmp_path / "output"

        orch = WithinSpeciesOrchestrator("test_sp", abundance_path, metadata_path, output_dir)
        orch.load_data()
        result = orch.run_differential_expression("nonexistent_column")
        assert result is None

    def test_run_all(self, tmp_path, synthetic_abundance, synthetic_metadata):
        """Test full run_all pipeline executes without error."""
        abundance_path, _ = synthetic_abundance
        metadata_path, _ = synthetic_metadata
        output_dir = tmp_path / "output"

        orch = WithinSpeciesOrchestrator("test_sp", abundance_path, metadata_path, output_dir)
        orch.run_all(condition_cols=["tissue", "caste", "sex", "nonexistent"])

        # PCA output should exist
        assert (output_dir / "test_sp_pca_coordinates.tsv").exists()
        # DE for tissue should exist (2 groups)
        assert (output_dir / "test_sp_DE_tissue.tsv").exists()
        # DE for caste should exist (2 groups in full dataset)
        assert (output_dir / "test_sp_DE_caste.tsv").exists()

    def test_load_data_no_common_samples_raises(self, tmp_path):
        """Test that mismatched samples raise ValueError."""
        # Create abundance with different sample names
        abundance = pd.DataFrame(
            np.random.poisson(50, (10, 3)),
            index=[f"g{i}" for i in range(10)],
            columns=["X1", "X2", "X3"],
        )
        meta = pd.DataFrame({"run": ["Y1", "Y2", "Y3"], "tissue": ["a", "b", "c"]})

        a_path = tmp_path / "ab.tsv"
        m_path = tmp_path / "me.tsv"
        abundance.to_csv(a_path, sep="\t")
        meta.to_csv(m_path, sep="\t", index=False)

        orch = WithinSpeciesOrchestrator("fail_sp", a_path, m_path, tmp_path / "out")
        with pytest.raises(ValueError, match="No common samples"):
            orch.load_data()

    def test_pca_before_load_raises(self, tmp_path, synthetic_abundance, synthetic_metadata):
        """Test that calling PCA before loading data raises ValueError."""
        abundance_path, _ = synthetic_abundance
        metadata_path, _ = synthetic_metadata
        output_dir = tmp_path / "output"

        orch = WithinSpeciesOrchestrator("test_sp", abundance_path, metadata_path, output_dir)
        with pytest.raises(ValueError, match="Data not loaded"):
            orch.run_pca()


# =============================================================================
# AcrossSpeciesOrchestrator Tests
# =============================================================================


class TestAcrossSpeciesOrchestrator:
    """Test suite for AcrossSpeciesOrchestrator."""

    def test_instantiation(self, tmp_path, ortholog_table):
        """Test that orchestrator can be instantiated."""
        orth_path, _ = ortholog_table
        orch = AcrossSpeciesOrchestrator(orth_path, tmp_path / "output")
        assert orch.output_dir.exists()

    def test_load_orthologs(self, tmp_path, ortholog_table):
        """Test ortholog table loading."""
        orth_path, _ = ortholog_table
        orch = AcrossSpeciesOrchestrator(orth_path, tmp_path / "output")
        orch.load_orthologs()

        assert orch.ortholog_df is not None
        assert len(orch.ortholog_df) == 50
        assert "Orthogroup" in orch.ortholog_df.columns

    def test_load_species_expression(self, tmp_path, ortholog_table, three_species_data):
        """Test species expression loading."""
        orth_path, _ = ortholog_table
        orch = AcrossSpeciesOrchestrator(orth_path, tmp_path / "output")

        for sp_name, (sp_path, _) in three_species_data.items():
            orch.load_species_expression(sp_name, sp_path)

        assert len(orch.species_expressions) == 3
        for sp_name in three_species_data:
            assert sp_name in orch.species_expressions

    def test_generate_pairwise_maps(self, tmp_path, ortholog_table, three_species_data):
        """Test pairwise ortholog map generation."""
        orth_path, _ = ortholog_table
        orch = AcrossSpeciesOrchestrator(orth_path, tmp_path / "output")
        orch.load_orthologs()

        for sp_name, (sp_path, _) in three_species_data.items():
            orch.load_species_expression(sp_name, sp_path)

        maps = orch.generate_pairwise_maps()

        # 3 species -> 6 directional pairs (A->B, B->A, etc.)
        assert len(maps) == 6
        for key, orth_map in maps.items():
            assert isinstance(key, tuple)
            assert len(key) == 2
            assert isinstance(orth_map, dict)
            assert len(orth_map) > 0

    def test_run_comparative_analysis(self, tmp_path, ortholog_table, three_species_data):
        """Test full comparative analysis produces output files."""
        orth_path, _ = ortholog_table
        output_dir = tmp_path / "output"
        orch = AcrossSpeciesOrchestrator(orth_path, output_dir)
        orch.load_orthologs()

        for sp_name, (sp_path, _) in three_species_data.items():
            orch.load_species_expression(sp_name, sp_path)

        orch.run_comparative_analysis()

        # Check output files were created
        assert (output_dir / "conservation_scores.tsv").exists()
        assert (output_dir / "expression_divergence_matrix.tsv").exists()

        # Validate conservation scores
        cons = pd.read_csv(output_dir / "conservation_scores.tsv", sep="\t")
        assert "gene_id" in cons.columns
        assert "mean_conservation" in cons.columns
        assert len(cons) > 0

        # Validate divergence matrix
        div = pd.read_csv(output_dir / "expression_divergence_matrix.tsv", sep="\t", index_col=0)
        assert div.shape[0] == div.shape[1]  # Square matrix
        assert div.shape[0] == 3  # 3 species
        # Diagonal should be ~0
        for sp in div.index:
            assert abs(div.loc[sp, sp]) < 1e-10

    def test_no_expression_data_raises(self, tmp_path, ortholog_table):
        """Test that running analysis without loading expression raises error."""
        orth_path, _ = ortholog_table
        orch = AcrossSpeciesOrchestrator(orth_path, tmp_path / "output")
        with pytest.raises(ValueError, match="No species expression"):
            orch.run_comparative_analysis()

    def test_conservation_scores_range(self, tmp_path, ortholog_table, three_species_data):
        """Test that conservation scores are within expected ranges."""
        orth_path, _ = ortholog_table
        output_dir = tmp_path / "output"
        orch = AcrossSpeciesOrchestrator(orth_path, output_dir)
        orch.load_orthologs()

        for sp_name, (sp_path, _) in three_species_data.items():
            orch.load_species_expression(sp_name, sp_path)

        orch.run_comparative_analysis()

        cons = pd.read_csv(output_dir / "conservation_scores.tsv", sep="\t")
        # mean_conservation should be between -1 and 1
        assert cons["mean_conservation"].min() >= -1.0
        assert cons["mean_conservation"].max() <= 1.0

    def test_divergence_matrix_symmetry(self, tmp_path, ortholog_table, three_species_data):
        """Test that the divergence matrix is symmetric."""
        orth_path, _ = ortholog_table
        output_dir = tmp_path / "output"
        orch = AcrossSpeciesOrchestrator(orth_path, output_dir)
        orch.load_orthologs()

        for sp_name, (sp_path, _) in three_species_data.items():
            orch.load_species_expression(sp_name, sp_path)

        orch.run_comparative_analysis()

        div = pd.read_csv(output_dir / "expression_divergence_matrix.tsv", sep="\t", index_col=0)
        # Check symmetry
        for sp_a in div.index:
            for sp_b in div.columns:
                assert abs(div.loc[sp_a, sp_b] - div.loc[sp_b, sp_a]) < 1e-10


# =============================================================================
# Cross-Module Integration: Within + Across
# =============================================================================


class TestWithinToAcrossIntegration:
    """Test that within-species outputs feed correctly into across-species analysis."""

    def test_pca_output_is_loadable_as_expression(self, tmp_path, synthetic_abundance, synthetic_metadata):
        """Test that within-species PCA coordinates can be read back as a data frame."""
        abundance_path, _ = synthetic_abundance
        metadata_path, _ = synthetic_metadata
        output_dir = tmp_path / "output"

        orch = WithinSpeciesOrchestrator("sp1", abundance_path, metadata_path, output_dir)
        orch.load_data()
        orch.run_pca()

        # Read back the PCA results
        pca_file = output_dir / "sp1_pca_coordinates.tsv"
        pca_df = pd.read_csv(pca_file, sep="\t", index_col=0)
        assert pca_df.shape[0] == 12  # samples
        assert pca_df.shape[1] == 2  # components

    def test_de_output_has_standard_columns(self, tmp_path, synthetic_abundance, synthetic_metadata):
        """Test that DE output TSV has standard columns for downstream processing."""
        abundance_path, _ = synthetic_abundance
        metadata_path, _ = synthetic_metadata
        output_dir = tmp_path / "output"

        orch = WithinSpeciesOrchestrator("sp1", abundance_path, metadata_path, output_dir)
        orch.load_data()
        result = orch.run_differential_expression("tissue")

        required_cols = {"log2_fold_change", "p_value", "regulation"}
        assert required_cols.issubset(set(result.columns))
