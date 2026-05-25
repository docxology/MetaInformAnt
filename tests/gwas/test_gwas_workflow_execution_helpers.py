"""Tests for private GWAS workflow execution helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.gwas.workflow.workflow_execution import (
    _apply_qvalue_and_bonferroni_annotations,
    _lambda_gc_from_associations,
    _metadata_path_from_config,
    _sample_major_to_variant_major,
    _validated_gwas_io_paths,
    _variant_major_to_sample_major,
    _write_summary_outputs,
)


def test_validated_gwas_io_paths_creates_output_dir(tmp_path: Path) -> None:
    """Input validation should use real files and create output directories."""
    vcf_path = tmp_path / "input.vcf"
    phenotype_path = tmp_path / "phenotypes.tsv"
    vcf_path.write_text("##fileformat=VCFv4.2\n")
    phenotype_path.write_text("sample\ttrait\nS1\t1.0\n")

    out_dir = tmp_path / "results"
    resolved_vcf, resolved_pheno, resolved_output = _validated_gwas_io_paths(vcf_path, phenotype_path, out_dir)

    assert resolved_vcf == vcf_path
    assert resolved_pheno == phenotype_path
    assert resolved_output == out_dir
    assert out_dir.is_dir()


def test_validated_gwas_io_paths_rejects_missing_inputs(tmp_path: Path) -> None:
    """Missing input files should fail before workflow execution starts."""
    phenotype_path = tmp_path / "phenotypes.tsv"
    phenotype_path.write_text("sample\ttrait\nS1\t1.0\n")

    with pytest.raises(FileNotFoundError):
        _validated_gwas_io_paths(tmp_path / "missing.vcf", phenotype_path, tmp_path / "results")


def test_genotype_matrix_transposes_round_trip() -> None:
    """GWAS helpers should preserve sample/variant orientation exactly."""
    sample_major = [[0, 1, 2], [2, 1, 0]]

    variant_major, n_samples, n_variants = _sample_major_to_variant_major(sample_major)

    assert variant_major == [[0, 2], [1, 1], [2, 0]]
    assert (n_samples, n_variants) == (2, 3)
    assert _variant_major_to_sample_major(variant_major, n_samples) == sample_major


def test_metadata_path_from_config_supports_flat_and_nested() -> None:
    """Flat metadata_file takes precedence, while nested samples metadata remains supported."""
    assert _metadata_path_from_config({"metadata_file": "flat.tsv"}) == "flat.tsv"
    assert _metadata_path_from_config({"samples": {"metadata_file": "nested.tsv"}}) == "nested.tsv"
    assert _metadata_path_from_config({"samples": []}) is None


def test_correction_helpers_annotate_association_rows() -> None:
    """Correction helper should mutate rows with q-values and Bonferroni flags."""
    rows = [
        {"p_value": 0.001},
        {"p_value": 0.2},
        {"p_value": 0.8},
    ]

    _apply_qvalue_and_bonferroni_annotations(rows)

    assert all("q_value" in row for row in rows)
    assert all("bonferroni_significant" in row for row in rows)
    assert _lambda_gc_from_associations(rows) is not None


def test_write_summary_outputs_uses_standard_files(tmp_path: Path) -> None:
    """Summary output helper should write all standard GWAS report artifacts."""
    assoc_results = [{"beta": 0.5, "se": 0.1, "p_value": 1e-6, "q_value": 2e-6, "n_samples": 10, "maf": 0.3}]
    variant_info = [{"chrom": "chr1", "pos": 100, "id": "rs1", "ref": "A", "alt": "G"}]

    outputs = _write_summary_outputs(assoc_results, variant_info, tmp_path, threshold=1e-5)

    assert Path(outputs["summary_stats_path"]).exists()
    assert Path(outputs["significant_hits_path"]).exists()
    assert (tmp_path / "results_summary.json").exists()
    assert outputs["summary"]["n_variants_tested"] == 1
