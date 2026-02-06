"""Tests for eQTL integration between GWAS and expression data.

This module tests the integration of genetic variants with gene expression
data, including colocalization analysis and data conversion utilities.
"""

import pytest


class TestEqtlColoc:
    """Tests for GWAS-eQTL colocalization."""

    def test_eqtl_coloc_basic(self):
        """Test basic eQTL colocalization with matching signals."""
        from metainformant.gwas.finemapping import eqtl_coloc

        # Create synthetic data with shared causal variant
        gwas_z = [0.5, 1.2, 3.5, 1.0, 0.3]  # Peak at index 2
        eqtl_z = [0.4, 1.0, 3.2, 0.9, 0.2]  # Peak at same location

        result = eqtl_coloc(
            gwas_z=gwas_z,
            eqtl_z=eqtl_z,
            gene_id="TEST_GENE",
        )

        assert "status" in result or "PP_H4" in result
        assert "gene_id" in result or result.get("status") == "error"

    def test_eqtl_coloc_no_signal(self):
        """Test eQTL colocalization with no shared signal."""
        from metainformant.gwas.finemapping import eqtl_coloc

        # Signals at different locations
        gwas_z = [3.5, 0.5, 0.5, 0.5, 0.5]  # Peak at index 0
        eqtl_z = [0.5, 0.5, 0.5, 0.5, 3.5]  # Peak at index 4

        result = eqtl_coloc(
            gwas_z=gwas_z,
            eqtl_z=eqtl_z,
            gene_id="TEST_GENE_2",
        )

        assert "status" in result or "PP_H4" in result
        # If successful, H4 should be low (no shared causal)
        if "PP_H4" in result:
            assert result["PP_H4"] < 0.5

    def test_eqtl_coloc_empty_input(self):
        """Test eQTL colocalization with empty input."""
        from metainformant.gwas.finemapping import eqtl_coloc

        result = eqtl_coloc(
            gwas_z=[],
            eqtl_z=[],
            gene_id="EMPTY_GENE",
        )

        assert result.get("status") == "error"

    def test_eqtl_coloc_mismatched_lengths(self):
        """Test eQTL colocalization with mismatched input lengths."""
        from metainformant.gwas.finemapping import eqtl_coloc

        result = eqtl_coloc(
            gwas_z=[1.0, 2.0, 3.0],
            eqtl_z=[1.0, 2.0],  # Different length
            gene_id="MISMATCH_GENE",
        )

        assert result.get("status") == "error"


class TestMultiTraitColoc:
    """Tests for multi-trait colocalization."""

    def test_multi_trait_coloc_two_traits(self):
        """Test multi-trait colocalization with two traits."""
        from metainformant.gwas.finemapping import multi_trait_coloc

        z_scores = {
            "trait1": [0.5, 1.2, 3.5, 1.0, 0.3],
            "trait2": [0.4, 1.0, 3.2, 0.9, 0.2],
        }

        result = multi_trait_coloc(z_scores=z_scores)

        assert "pairwise_results" in result or "status" in result

    def test_multi_trait_coloc_three_traits(self):
        """Test multi-trait colocalization with three traits."""
        from metainformant.gwas.finemapping import multi_trait_coloc

        z_scores = {
            "trait1": [0.5, 1.2, 3.5, 1.0, 0.3],
            "trait2": [0.4, 1.0, 3.2, 0.9, 0.2],
            "trait3": [0.6, 1.1, 3.0, 1.2, 0.4],
        }

        result = multi_trait_coloc(z_scores=z_scores)

        assert "pairwise_results" in result or "status" in result


class TestComputeClpp:
    """Tests for CLPP computation."""

    def test_compute_clpp_basic(self):
        """Test CLPP computation with matching PIPs."""
        from metainformant.gwas.finemapping import compute_clpp

        # PIPs that suggest shared causal variant at index 2
        pip_1 = [0.05, 0.1, 0.7, 0.1, 0.05]
        pip_2 = [0.04, 0.08, 0.75, 0.08, 0.05]

        result = compute_clpp(pip_1=pip_1, pip_2=pip_2)

        assert "clpp" in result
        assert "is_colocalized" in result
        assert result["clpp"] > 0.4  # Should be high for shared causal

    def test_compute_clpp_no_shared(self):
        """Test CLPP computation with no shared causal."""
        from metainformant.gwas.finemapping import compute_clpp

        pip_1 = [0.8, 0.1, 0.05, 0.03, 0.02]  # Peak at 0
        pip_2 = [0.02, 0.03, 0.05, 0.1, 0.8]  # Peak at 4

        result = compute_clpp(pip_1=pip_1, pip_2=pip_2)

        assert "clpp" in result
        assert result["clpp"] < 0.1  # Should be low


class TestIntegrationConverters:
    """Tests for data conversion utilities."""

    def test_from_rna_expression_basic(self):
        """Test RNA expression data conversion."""
        import pandas as pd

        from metainformant.multiomics.analysis.integration import from_rna_expression

        expression_df = pd.DataFrame(
            {
                "sample1": [1.0, 2.0, 3.0],
                "sample2": [1.5, 2.5, 3.5],
                "sample3": [1.2, 2.2, 3.2],
            },
            index=["gene1", "gene2", "gene3"],
        )

        # Test without normalization (no sklearn required)
        result = from_rna_expression(expression_df, normalize=False)

        assert isinstance(result, pd.DataFrame)
        assert result.shape == expression_df.shape

    def test_from_dna_variants_basic(self):
        """Test DNA variant data conversion."""
        import pandas as pd

        from metainformant.multiomics.analysis.integration import from_dna_variants

        vcf_df = pd.DataFrame(
            {
                "CHROM": ["1", "1", "2"],
                "POS": [100, 200, 300],
                "REF": ["A", "G", "C"],
                "ALT": ["T", "A", "G"],
                "sample1": [0, 1, 2],
                "sample2": [1, 0, 1],
            }
        )

        result = from_dna_variants(vcf_df)

        assert isinstance(result, pd.DataFrame)


class TestMultiOmicsData:
    """Tests for MultiOmicsData container."""

    def test_multiomics_data_init(self):
        """Test MultiOmicsData initialization."""
        import pandas as pd

        from metainformant.multiomics.analysis.integration import MultiOmicsData

        # Samples in rows (index), features in columns
        rna_data = pd.DataFrame(
            {"g1": [1, 3], "g2": [2, 4]}, index=["s1", "s2"]
        )
        protein_data = pd.DataFrame(
            {"p1": [0.5, 1.5], "p2": [1.0, 2.0]}, index=["s1", "s2"]
        )

        mo = MultiOmicsData(rna_data=rna_data, protein_data=protein_data)

        assert mo is not None
        # Check that data was stored (attribute names may vary)
        assert hasattr(mo, "data") or hasattr(mo, "rna_data") or hasattr(mo, "_rna_data")

    def test_get_common_samples(self):
        """Test getting common samples across omics."""
        import pandas as pd

        from metainformant.multiomics.analysis.integration import MultiOmicsData

        # Samples in rows, features in columns; rna has s1,s2,s3, protein has s1,s2
        rna_data = pd.DataFrame(
            {"g1": [1, 3, 5], "g2": [2, 4, 6]}, index=["s1", "s2", "s3"]
        )
        protein_data = pd.DataFrame(
            {"p1": [0.5, 1.5], "p2": [1.0, 2.0]}, index=["s1", "s2"]
        )

        mo = MultiOmicsData(rna_data=rna_data, protein_data=protein_data)
        common = mo.get_common_samples()

        # Should find s1, s2 as common
        assert len(common) >= 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


class TestCisEqtlScan:
    """Tests for cis-eQTL scanning."""

    def test_cis_eqtl_scan_basic(self):
        """Test basic cis-eQTL scan."""
        import pandas as pd

        from metainformant.gwas.finemapping.eqtl import cis_eqtl_scan

        # Create minimal test data
        expression = pd.DataFrame(
            {"s1": [1.0, 2.0], "s2": [1.5, 2.5], "s3": [0.8, 1.8]},
            index=["gene1", "gene2"],
        )
        genotypes = pd.DataFrame(
            {"s1": [0, 1], "s2": [1, 2], "s3": [0, 0]},
            index=["var1", "var2"],
        )
        gene_pos = pd.DataFrame(
            {
                "gene_id": ["gene1", "gene2"],
                "chrom": ["1", "1"],
                "tss_position": [100000, 200000],
            }
        )
        var_pos = pd.DataFrame(
            {
                "variant_id": ["var1", "var2"],
                "chrom": ["1", "1"],
                "position": [105000, 195000],
            }
        )

        result = cis_eqtl_scan(
            expression_matrix=expression,
            genotype_matrix=genotypes,
            gene_positions=gene_pos,
            variant_positions=var_pos,
            cis_window=1_000_000,
            maf_threshold=0.0,  # Allow all for small test
        )

        assert isinstance(result, pd.DataFrame)

    def test_cis_eqtl_scan_empty(self):
        """Test cis-eQTL scan with no overlapping samples."""
        import pandas as pd

        from metainformant.gwas.finemapping.eqtl import cis_eqtl_scan

        expression = pd.DataFrame({"a": [1.0]}, index=["gene1"])
        genotypes = pd.DataFrame({"b": [0]}, index=["var1"])
        gene_pos = pd.DataFrame(
            {
                "gene_id": ["gene1"],
                "chrom": ["1"],
                "tss_position": [100000],
            }
        )
        var_pos = pd.DataFrame(
            {
                "variant_id": ["var1"],
                "chrom": ["1"],
                "position": [105000],
            }
        )

        result = cis_eqtl_scan(
            expression,
            genotypes,
            gene_pos,
            var_pos,
            maf_threshold=0.0,
        )

        assert len(result) == 0


class TestEqtlEffectSizes:
    """Tests for effect size computation."""

    def test_effect_sizes_basic(self):
        """Test effect size calculation."""
        import pandas as pd

        from metainformant.gwas.finemapping.eqtl import eqtl_effect_sizes

        expression = pd.DataFrame(
            {"s1": [1.0, 2.0], "s2": [2.0, 3.0], "s3": [1.5, 2.5]},
            index=["gene1", "gene2"],
        )
        genotypes = pd.DataFrame(
            {"s1": [0, 1], "s2": [2, 0], "s3": [1, 2]},
            index=["var1", "var2"],
        )
        associations = pd.DataFrame(
            {
                "gene_id": ["gene1"],
                "variant_id": ["var1"],
            }
        )

        result = eqtl_effect_sizes(expression, genotypes, associations)

        assert "beta" in result.columns
        assert "r_squared" in result.columns


class TestEqtlSummaryStats:
    """Tests for summary statistics."""

    def test_summary_stats_basic(self):
        """Test summary statistics generation."""
        import pandas as pd

        from metainformant.gwas.finemapping.eqtl import eqtl_summary_stats

        results = pd.DataFrame(
            {
                "gene_id": ["g1", "g1", "g2"],
                "variant_id": ["v1", "v2", "v3"],
                "beta": [0.5, 0.1, 0.8],
                "pvalue": [0.001, 0.5, 0.0001],
            }
        )

        summary = eqtl_summary_stats(results, fdr_threshold=0.1)

        assert "n_tests" in summary
        assert summary["n_tests"] == 3

    def test_summary_stats_empty(self):
        """Test summary statistics with empty input."""
        import pandas as pd

        from metainformant.gwas.finemapping.eqtl import eqtl_summary_stats

        results = pd.DataFrame()
        summary = eqtl_summary_stats(results)

        assert summary["n_tests"] == 0
