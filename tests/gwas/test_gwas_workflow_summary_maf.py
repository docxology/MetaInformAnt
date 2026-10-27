"""Regression tests: the GWAS workflows write the MAF column of their summary statistics.

Defect: ``write_summary_statistics`` / ``write_significant_hits`` take an optional ``genotypes=`` matrix
(variants x samples, 0/1/2 dosages) from which they compute the MAF of every row. The workflow call sites in
``workflow_execution`` never passed it, so ``summary_statistics.tsv`` and ``significant_hits.tsv`` of every pipeline
run carried ``NA`` in the MAF column. (The writers themselves are covered by
``test_gwas_summary_stats_missing_values.py``; what was missing is the hand-over from the pipelines.)

Oracles: the VCF is written by hand and the expected minor allele frequencies are counted from its genotype text
with plain string handling, independent of ``parse_vcf_full`` and of the writers' MAF code. The same values are also
worked out by hand in ``HAND_COUNTED_MAF``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from metainformant.gwas.workflow.workflow_execution import (
    _write_summary_outputs,
    execute_gwas_workflow,
    run_gwas,
    run_multi_trait_gwas,
)

# Eight samples, six variants. rsE has one missing call; rsF is heterozygous in every sample.
VCF_TEXT = (
    "##fileformat=VCFv4.2\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\n"
    "chr1\t100\trsA\tA\tG\t60\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/1\t1/1\t0/0\t0/1\t0/0\n"
    "chr1\t200\trsB\tC\tT\t60\tPASS\t.\tGT\t0/1\t0/1\t1/1\t0/1\t0/0\t1/1\t0/1\t1/1\n"
    "chr1\t300\trsC\tG\tA\t60\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0\t0/1\n"
    "chr1\t400\trsD\tT\tC\t60\tPASS\t.\tGT\t1/1\t0/1\t0/1\t0/0\t1/1\t0/1\t1/1\t0/1\n"
    "chr2\t100\trsE\tA\tC\t60\tPASS\t.\tGT\t0/1\t./.\t0/0\t0/1\t0/0\t0/0\t1/1\t0/1\n"
    "chr2\t200\trsF\tG\tT\t60\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\n"
)
PHENOTYPES_TEXT = (
    "sample_id\tyield\tvigor\n"
    "S1\t1.2\t3.1\nS2\t2.3\t2.9\nS3\t2.9\t4.4\nS4\t1.8\t3.3\nS5\t3.6\t5.0\nS6\t2.1\t3.9\nS7\t3.9\t4.1\nS8\t1.1\t2.2\n"
)
# The same trait for the first six samples only: the workflows then analyse the first six samples.
PHENOTYPES_FIRST_SIX_TEXT = "sample_id\tyield\nS1\t1.2\nS2\t2.3\nS3\t2.9\nS4\t1.8\nS5\t3.6\nS6\t2.1\n"

# Pass-through quality control: no variant is removed, so every VCF line becomes a row of the output.
KEEP_ALL = {"min_maf": 0.0, "max_missing": 1.0, "min_hwe_p": 0.0}

# Alternate alleles counted by hand over the eight samples.
#   rsA 0+1+0+1+2+0+1+0 = 5 of 16 -> 0.3125     rsB 1+1+2+1+0+2+1+2 = 10 of 16 -> alt 0.625 -> MAF 0.375
#   rsC 0+0+0+0+1+0+0+1 = 2 of 16 -> 0.125      rsD 2+1+1+0+2+1+2+1 = 10 of 16 -> MAF 0.375
#   rsE seven called samples: 1+0+1+0+0+2+1 = 5 of 14 -> 0.357142...   rsF 8 of 16 -> 0.5
HAND_COUNTED_MAF = {"rsA": 0.3125, "rsB": 0.375, "rsC": 0.125, "rsD": 0.375, "rsE": 5 / 14, "rsF": 0.5}


def maf_by_snp_from_vcf_text(n_samples: int | None = None) -> dict[str, float]:
    """Minor allele frequency per SNP, counted from the genotype text of the first ``n_samples`` samples."""
    result: dict[str, float] = {}
    for line in VCF_TEXT.splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        calls = columns[9:] if n_samples is None else columns[9 : 9 + n_samples]
        alleles = [allele for call in calls for allele in call.replace("|", "/").split("/") if allele != "."]
        alt_frequency = sum(1 for allele in alleles if allele == "1") / len(alleles)
        result[columns[2]] = min(alt_frequency, 1.0 - alt_frequency)
    return result


def read_table(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().strip().split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def assert_maf_column(path: Path, expected: dict[str, float], *, minimum_rows: int) -> None:
    """Every row's MAF is the value counted from the genotypes (written with four decimals), never NA."""
    rows = read_table(path)
    assert len(rows) >= minimum_rows
    for row in rows:
        assert row["MAF"] != "NA", f"{row['SNP']} has no MAF in {path.name}"
        assert float(row["MAF"]) == pytest.approx(expected[row["SNP"]], abs=5.1e-5)


@pytest.fixture
def dataset(tmp_path: Path) -> dict[str, Path]:
    vcf = tmp_path / "mini.vcf"
    vcf.write_text(VCF_TEXT)
    phenotypes = tmp_path / "phenotypes.tsv"
    phenotypes.write_text(PHENOTYPES_TEXT)
    return {"vcf": vcf, "phenotypes": phenotypes}


def test_the_independent_oracle_agrees_with_the_hand_count() -> None:
    """Guard the oracle itself before it is used against the pipelines."""
    counted = maf_by_snp_from_vcf_text()

    assert counted.keys() == HAND_COUNTED_MAF.keys()
    for snp, expected in HAND_COUNTED_MAF.items():
        assert counted[snp] == pytest.approx(expected, abs=1e-12)


def test_run_gwas_writes_the_maf_of_every_variant(tmp_path: Path, dataset: dict[str, Path]) -> None:
    out = tmp_path / "run_gwas"

    result = run_gwas(
        dataset["vcf"],
        dataset["phenotypes"],
        {"model": "linear", "trait": "yield", "quality_control": KEEP_ALL, "significance_threshold": 1.0},
        output_dir=out,
    )

    assert result["status"] == "success"
    expected = maf_by_snp_from_vcf_text()
    assert_maf_column(out / "summary_statistics.tsv", expected, minimum_rows=6)
    # significance_threshold 1.0 puts every test with p < 1 into the hits table, which needs the MAF too.
    assert_maf_column(out / "significant_hits.tsv", expected, minimum_rows=3)


def test_run_gwas_maf_counts_the_samples_the_tests_used(tmp_path: Path, dataset: dict[str, Path]) -> None:
    """With phenotypes for six of the eight samples the tests run on the first six; so must the MAF."""
    phenotypes = tmp_path / "phenotypes_first_six.tsv"
    phenotypes.write_text(PHENOTYPES_FIRST_SIX_TEXT)
    out = tmp_path / "run_gwas_six"

    result = run_gwas(
        dataset["vcf"],
        phenotypes,
        {"model": "linear", "trait": "yield", "quality_control": KEEP_ALL, "significance_threshold": 1.0},
        output_dir=out,
    )

    assert result["status"] == "success"
    expected_six = maf_by_snp_from_vcf_text(n_samples=6)
    assert expected_six != maf_by_snp_from_vcf_text()  # the two sample sets give different MAFs
    assert_maf_column(out / "summary_statistics.tsv", expected_six, minimum_rows=6)
    assert {row["N"] for row in read_table(out / "summary_statistics.tsv")} <= {"5", "6"}  # rsE misses one call
    assert_maf_column(out / "significant_hits.tsv", expected_six, minimum_rows=3)


def test_execute_gwas_workflow_writes_the_maf_of_every_retained_variant(
    tmp_path: Path, dataset: dict[str, Path]
) -> None:
    """The flat-config workflow LD-prunes the variants; the genotype rows must stay aligned with the results."""
    out = tmp_path / "execute"
    config: dict[str, Any] = {
        "vcf_path": str(dataset["vcf"]),
        "phenotype_path": str(dataset["phenotypes"]),
        "output_dir": str(out),
        "results_dir": str(out),
        "trait": "yield",
        "model": "linear",
        "quality_control": KEEP_ALL,
        "significance_threshold": 1.0,
    }

    result = execute_gwas_workflow(config)

    assert result["status"] == "completed"
    assert "summary_statistics" in result["steps_completed"]
    rows = read_table(out / "summary_statistics.tsv")
    assert 2 <= len(rows) <= 6
    assert_maf_column(out / "summary_statistics.tsv", maf_by_snp_from_vcf_text(), minimum_rows=2)
    # The hits file is written from the same call site and must carry the column as well (it may have no rows).
    for row in read_table(out / "significant_hits.tsv"):
        assert row["MAF"] != "NA"


def test_run_multi_trait_gwas_writes_the_maf_for_every_trait(tmp_path: Path, dataset: dict[str, Path]) -> None:
    out = tmp_path / "multi"

    result = run_multi_trait_gwas(
        dataset["vcf"], dataset["phenotypes"], ["yield", "vigor"], {"model": "linear", "quality_control": KEEP_ALL}, out
    )

    assert result["status"] == "success"
    expected = maf_by_snp_from_vcf_text()
    for trait in ("yield", "vigor"):
        assert result["trait_results"][trait]["status"] == "success"
        assert_maf_column(out / trait / "summary_statistics.tsv", expected, minimum_rows=6)


def test_write_summary_outputs_passes_the_genotypes_to_both_tables(tmp_path: Path) -> None:
    """The helper shared by ``run_gwas``: genotypes (variants x samples) fill the MAF of both output tables."""
    results = [
        {"beta": 0.5, "se": 0.1, "p_value": 1e-9, "n_samples": 4},
        {"beta": 0.1, "se": 0.1, "p_value": 0.5, "n_samples": 4},
    ]
    variants = [
        {"chrom": "1", "pos": 10, "id": "rs1", "ref": "A", "alt": ["G"]},
        {"chrom": "1", "pos": 20, "id": "rs2", "ref": "C", "alt": ["T"]},
    ]
    genotypes = [[0, 0, 1, 2], [2, 2, 2, 1]]  # rs1: 3 of 8 -> 0.375; rs2: 7 of 8 -> alt 0.875 -> MAF 0.125

    _write_summary_outputs(results, variants, tmp_path, threshold=1e-6, genotypes=genotypes)

    assert [row["MAF"] for row in read_table(tmp_path / "summary_statistics.tsv")] == ["0.3750", "0.1250"]
    assert [row["MAF"] for row in read_table(tmp_path / "significant_hits.tsv")] == ["0.3750"]


def test_write_summary_outputs_without_genotypes_still_reports_na(tmp_path: Path) -> None:
    """Nothing is invented when the caller has no genotypes: the MAF stays NA."""
    results = [{"beta": 0.5, "se": 0.1, "p_value": 1e-9, "n_samples": 4}]
    variants = [{"chrom": "1", "pos": 10, "id": "rs1", "ref": "A", "alt": ["G"]}]

    _write_summary_outputs(results, variants, tmp_path, threshold=1e-6)

    assert [row["MAF"] for row in read_table(tmp_path / "summary_statistics.tsv")] == ["NA"]
