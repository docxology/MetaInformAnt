"""Tests for RNA species discovery and configuration generation functions.

This module tests discovery functionality following NO_MOCKING_POLICY.
Network-dependent tests are skipped when offline.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.rna import discovery


class TestGenerateConfigYaml:
    """Test generate_config_yaml function."""

    def test_generate_config_yaml_basic(self):
        """Test generate_config_yaml generates basic config."""
        species_name = "Homo sapiens"
        species_data = {
            "sample_count": 100,
            "taxonomy_id": "9606",
            "run_ids": ["SRR123", "SRR456"],
            "study_ids": ["SRP001"],
        }
        
        yaml_content = discovery.generate_config_yaml(species_name, species_data)
        assert isinstance(yaml_content, str)
        assert "Homo sapiens" in yaml_content
        assert "9606" in yaml_content
        assert "work_dir" in yaml_content
        assert "species_list" in yaml_content

    def test_generate_config_yaml_with_genome_info(self):
        """Test generate_config_yaml with genome information."""
        species_name = "Homo sapiens"
        species_data = {
            "sample_count": 100,
            "taxonomy_id": "9606",
            "run_ids": ["SRR123"],
            "study_ids": ["SRP001"],
        }
        genome_info = {
            "accession": "GCF_000001405.40",
            "assembly_name": "GRCh38.p14",
            "level": "chromosome",
            "annotation_release": "110",
        }
        
        yaml_content = discovery.generate_config_yaml(species_name, species_data, genome_info)
        assert "GCF_000001405.40" in yaml_content
        assert "genome:" in yaml_content
        assert "ftp_url" in yaml_content

    def test_generate_config_yaml_with_repo_root(self):
        """Test generate_config_yaml with repo_root parameter."""
        species_name = "Mus musculus"
        species_data = {
            "sample_count": 50,
            "taxonomy_id": "10090",
            "run_ids": ["SRR789"],
            "study_ids": ["SRP002"],
        }
        repo_root = Path("/test/repo")
        
        yaml_content = discovery.generate_config_yaml(species_name, species_data, repo_root=repo_root)
        assert "/test/repo" in yaml_content
        assert "output/amalgkit" in yaml_content


class TestSelectBestAssembly:
    """Test _select_best_assembly function."""

    def test_select_best_assembly_empty(self):
        """Test _select_best_assembly with empty list."""
        result = discovery._select_best_assembly([])
        assert result is None

    def test_select_best_assembly_prefers_refseq(self):
        """Test _select_best_assembly prefers RefSeq over GenBank."""
        assemblies = [
            {
                "assembly": {
                    "assembly_accession": "GCA_000001405.1",
                    "assembly_level": "chromosome",
                },
                "assembly_stats": {"contig_n50": 1000000},
            },
            {
                "assembly": {
                    "assembly_accession": "GCF_000001405.1",
                    "assembly_level": "chromosome",
                },
                "assembly_stats": {"contig_n50": 1000000},
            },
        ]
        
        result = discovery._select_best_assembly(assemblies)
        assert result is not None
        assert result["assembly"]["assembly_accession"].startswith("GCF_")

    def test_select_best_assembly_prefers_chromosome_level(self):
        """Test _select_best_assembly prefers chromosome-level assemblies."""
        assemblies = [
            {
                "assembly": {
                    "assembly_accession": "GCF_000001405.1",
                    "assembly_level": "contig",
                },
                "assembly_stats": {"contig_n50": 1000000},
            },
            {
                "assembly": {
                    "assembly_accession": "GCF_000001405.2",
                    "assembly_level": "chromosome",
                },
                "assembly_stats": {"contig_n50": 1000000},
            },
        ]
        
        result = discovery._select_best_assembly(assemblies)
        assert result is not None
        assert result["assembly"]["assembly_level"] == "chromosome"


class TestDiscoveryDocumentation:
    """Test that discovery functions have proper documentation."""

    def test_all_functions_have_docstrings(self):
        """Verify all discovery functions have docstrings."""
        functions = [
            discovery.generate_config_yaml,
            discovery._select_best_assembly,
        ]
        
        for func in functions:
            assert func.__doc__ is not None, f"{func.__name__} missing docstring"
            assert len(func.__doc__.strip()) > 0

    @pytest.mark.network
    def test_search_species_with_rnaseq_requires_biopython(self):
        """Test that search_species_with_rnaseq requires Biopython."""
        try:
            from metainformant.rna.discovery import BIOPYTHON_AVAILABLE
            if not BIOPYTHON_AVAILABLE:
                pytest.skip("Biopython not available")
            
            # If available, test basic functionality
            result = discovery.search_species_with_rnaseq(
                '"Homo sapiens"[Organism] AND RNA-Seq[Strategy]',
                max_records=10,
            )
            assert isinstance(result, dict)
        except ImportError:
            pytest.skip("Biopython not available for SRA search")

    @pytest.mark.network
    def test_get_genome_info_requires_ncbi_datasets(self):
        """Test that get_genome_info requires ncbi-datasets-pylib."""
        try:
            from metainformant.rna.discovery import NCBI_DATASETS_AVAILABLE
            if not NCBI_DATASETS_AVAILABLE:
                pytest.skip("ncbi-datasets-pylib not available")
            
            # If available, test basic functionality
            result = discovery.get_genome_info("9606", "Homo sapiens")
            # May return None if no assemblies found, or dict if found
            assert result is None or isinstance(result, dict)
        except ImportError:
            pytest.skip("ncbi-datasets-pylib not available for genome info")

