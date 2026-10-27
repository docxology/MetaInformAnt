"""Regression tests: a square genotype matrix is not silently mis-oriented.

Defect: ``_extract_genotype_matrix`` transposed the sample-major genotypes of ``parse_vcf_full`` only when their shape
told them apart from variant-major ones (``n_samples != n_variants``). With as many variants as samples the matrix was
left untransposed, so ``execute_gwas_workflow`` ran every association test (and wrote every MAF) on a *sample's*
genotypes across variants while labelling it with a *variant's* id. Nothing failed: the output was simply wrong.

A square matrix cannot be oriented by its shape, so the shape is no longer asked to. The caller that knows the
layout says so (``layout="sample_major"``, used by ``execute_gwas_workflow``); without ``layout`` a square matrix is
still returned unchanged (existing behaviour, now documented) but a warning is logged.

Oracles: the genotype text of a hand-written VCF is counted with plain string handling; the transposition of the small
matrices is written out by hand.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pytest

from metainformant.gwas.workflow.workflow_config import _extract_genotype_matrix
from metainformant.gwas.workflow.workflow_execution import execute_gwas_workflow, run_gwas

# Six samples and six variants: a square genotype matrix. Alternate alleles per variant: 4, 5, 2, 3, 1 and 6 of 12,
# so every variant has its own minor allele frequency (the frequencies per sample are different numbers).
SQUARE_VCF = (
    "##fileformat=VCFv4.2\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\n"
    "chr1\t100\trsA\tA\tG\t60\tPASS\t.\tGT\t0/1\t0/1\t1/1\t0/0\t0/0\t0/0\n"
    "chr1\t200\trsB\tC\tT\t60\tPASS\t.\tGT\t1/1\t0/1\t0/0\t0/1\t0/0\t0/1\n"
    "chr1\t300\trsC\tG\tA\t60\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/0\t0/1\t0/0\n"
    "chr2\t100\trsD\tT\tC\t60\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1\n"
    "chr2\t200\trsE\tA\tC\t60\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/1\n"
    "chr2\t300\trsF\tG\tT\t60\tPASS\t.\tGT\t1/1\t0/0\t1/1\t0/0\t0/1\t0/1\n"
)
PHENOTYPES = "sample_id\tyield\nS1\t1.2\nS2\t2.3\nS3\t2.9\nS4\t1.8\nS5\t3.6\nS6\t2.1\n"
KEEP_ALL = {"min_maf": 0.0, "max_missing": 1.0, "min_hwe_p": 0.0}
# 4/12, 5/12 (alt frequency 5/12), 2/12, 3/12, 1/12 and 6/12
HAND_COUNTED_MAF = {"rsA": 4 / 12, "rsB": 5 / 12, "rsC": 2 / 12, "rsD": 3 / 12, "rsE": 1 / 12, "rsF": 6 / 12}


def maf_by_snp(vcf_text: str) -> Dict[str, float]:
    """Minor allele frequency per SNP from the genotype text, independent of every parser of the package."""
    frequencies = {}
    for line in vcf_text.splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        alleles = [a for call in columns[9:] for a in call.split("/") if a != "."]
        alt = sum(1 for a in alleles if a == "1") / len(alleles)
        frequencies[columns[2]] = min(alt, 1 - alt)
    return frequencies


def read_table(path: Path) -> List[Dict[str, str]]:
    lines = path.read_text().strip().split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def vcf_record(samples: int, genotypes: Any, variants: int = 0) -> Dict[str, Any]:
    record: Dict[str, Any] = {"samples": [f"S{i}" for i in range(samples)], "genotypes": genotypes}
    if variants:
        record["variants"] = [{"id": f"v{i}"} for i in range(variants)]
    return record


def extract(record: Dict[str, Any], layout: Optional[str] = None) -> Any:
    """Call the helper the way callers did before the ``layout`` option existed when no layout is declared."""
    if layout is None:
        return _extract_genotype_matrix(record)
    return _extract_genotype_matrix(record, layout=layout)


# ---------------------------------------------------------------------------
# The helper
# ---------------------------------------------------------------------------

# Sample-major: row s holds the genotypes of sample s. Its transpose is variant-major.
SAMPLE_MAJOR_SQUARE = [[0, 1, 2], [0, 1, 2], [0, 1, 2]]
VARIANT_MAJOR_OF_THAT = [[0, 0, 0], [1, 1, 1], [2, 2, 2]]


def test_the_oracle_matrices_are_not_symmetric() -> None:
    """Guard the fixtures: a symmetric matrix could not tell the two orientations apart."""
    assert SAMPLE_MAJOR_SQUARE != VARIANT_MAJOR_OF_THAT
    assert np.array(SAMPLE_MAJOR_SQUARE).T.tolist() == VARIANT_MAJOR_OF_THAT


def test_a_declared_sample_major_square_matrix_is_transposed() -> None:
    result = _extract_genotype_matrix(vcf_record(3, SAMPLE_MAJOR_SQUARE, variants=3), layout="sample_major")

    assert result == VARIANT_MAJOR_OF_THAT


def test_a_declared_sample_major_matrix_is_transposed_without_variant_records_too() -> None:
    result = _extract_genotype_matrix(vcf_record(3, SAMPLE_MAJOR_SQUARE), layout="sample_major")

    assert result == VARIANT_MAJOR_OF_THAT


def test_a_declared_variant_major_square_matrix_is_returned_as_it_is() -> None:
    result = _extract_genotype_matrix(vcf_record(3, VARIANT_MAJOR_OF_THAT, variants=3), layout="variant_major")

    assert result == VARIANT_MAJOR_OF_THAT


def test_numpy_matrices_are_handled_like_lists() -> None:
    sample_major = np.array(SAMPLE_MAJOR_SQUARE)

    transposed = _extract_genotype_matrix(vcf_record(3, sample_major, variants=3), layout="sample_major")
    unchanged = _extract_genotype_matrix(vcf_record(3, sample_major.T, variants=3), layout="variant_major")

    assert isinstance(transposed, np.ndarray)
    assert transposed.tolist() == unchanged.tolist() == VARIANT_MAJOR_OF_THAT


@pytest.mark.parametrize("layout", ["sample_major", None])
def test_rectangular_sample_major_matrices_are_transposed_with_or_without_a_declaration(layout: Any) -> None:
    sample_major = [[0, 1, 2, 0], [1, 1, 0, 2], [2, 0, 1, 1]]  # 3 samples x 4 variants
    expected = [[0, 1, 2], [1, 1, 0], [2, 0, 1], [0, 2, 1]]

    assert extract(vcf_record(3, sample_major, variants=4), layout) == expected


@pytest.mark.parametrize("layout", ["variant_major", None])
def test_rectangular_variant_major_matrices_are_returned_as_they_are(layout: Any) -> None:
    variant_major = [[0, 1, 2], [1, 1, 0], [2, 0, 1], [0, 2, 1]]  # 4 variants x 3 samples

    assert extract(vcf_record(3, variant_major, variants=4), layout) == variant_major


def test_a_declaration_that_the_shape_contradicts_is_an_error() -> None:
    variant_major = [[0, 1, 2], [1, 1, 0], [2, 0, 1], [0, 2, 1]]  # 4 variants x 3 samples

    with pytest.raises(ValueError, match="layout='sample_major' needs one row per sample"):
        _extract_genotype_matrix(vcf_record(3, variant_major, variants=4), layout="sample_major")
    with pytest.raises(ValueError, match="layout='variant_major' needs one row per variant record"):
        _extract_genotype_matrix(vcf_record(3, [[0, 1, 2], [1, 1, 0]], variants=4), layout="variant_major")


def test_ragged_rows_are_reported_not_turned_into_an_index_error() -> None:
    with pytest.raises(ValueError, match="different lengths"):
        _extract_genotype_matrix(vcf_record(3, [[0, 1, 2], [1, 1]], variants=2), layout="sample_major")


def test_an_unknown_layout_is_rejected() -> None:
    with pytest.raises(ValueError, match="Unknown genotype layout 'rows'"):
        _extract_genotype_matrix(vcf_record(3, SAMPLE_MAJOR_SQUARE), layout="rows")


def test_without_a_declaration_a_square_matrix_is_returned_unchanged_with_a_warning(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """The documented assumption (variant-major) is kept, so existing callers behave as before, but not silently."""
    with caplog.at_level(logging.WARNING):
        result = _extract_genotype_matrix(vcf_record(3, SAMPLE_MAJOR_SQUARE, variants=3))

    assert result == SAMPLE_MAJOR_SQUARE
    assert "square" in caplog.text and "layout=" in caplog.text


def test_a_rectangular_matrix_gives_no_square_warning(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.WARNING):
        _extract_genotype_matrix(vcf_record(3, [[0, 1, 2, 0], [1, 1, 0, 2], [2, 0, 1, 1]], variants=4))

    assert "square" not in caplog.text


def test_the_docstring_states_the_assumption_for_square_matrices() -> None:
    doc = " ".join((_extract_genotype_matrix.__doc__ or "").split())

    assert "square matrix" in doc
    assert "assumed to be variant-major" in doc
    assert "layout" in doc


# ---------------------------------------------------------------------------
# The workflow that used to mislabel the data
# ---------------------------------------------------------------------------


@pytest.fixture
def square_dataset(tmp_path: Path) -> Dict[str, Path]:
    vcf = tmp_path / "square.vcf"
    vcf.write_text(SQUARE_VCF)
    phenotypes = tmp_path / "phenotypes.tsv"
    phenotypes.write_text(PHENOTYPES)
    return {"vcf": vcf, "phenotypes": phenotypes}


def test_the_hand_counted_frequencies_match_the_text_oracle() -> None:
    counted = maf_by_snp(SQUARE_VCF)

    for snp, expected in HAND_COUNTED_MAF.items():
        assert counted[snp] == pytest.approx(expected)


def test_execute_gwas_workflow_labels_each_variant_with_its_own_genotypes(
    tmp_path: Path, square_dataset: Dict[str, Path]
) -> None:
    """Six samples and six variants: every summary row must carry the MAF of its own SNP."""
    out = tmp_path / "execute"
    config = {
        "vcf_path": str(square_dataset["vcf"]),
        "phenotype_path": str(square_dataset["phenotypes"]),
        "output_dir": str(out),
        "results_dir": str(out),
        "trait": "yield",
        "model": "linear",
        "quality_control": KEEP_ALL,
        "ld_pruning": {"r2_threshold": 1.0},  # nothing is correlated beyond r2 = 1: keep every variant
        "structure": {"pca_components": 1},  # six samples cannot carry four PC covariates plus an intercept
    }

    result = execute_gwas_workflow(config)

    assert result["status"] == "completed"
    rows = read_table(out / "summary_statistics.tsv")
    assert [row["SNP"] for row in rows] == list(HAND_COUNTED_MAF)  # all six, in file order
    for row in rows:
        assert float(row["MAF"]) == pytest.approx(HAND_COUNTED_MAF[row["SNP"]], abs=5.1e-5), row["SNP"]
        assert row["N"] == "6"


def test_run_gwas_was_never_affected_and_still_agrees(tmp_path: Path, square_dataset: Dict[str, Path]) -> None:
    """``run_gwas`` always transposed (no shape heuristic): the same square data gives the same, correct MAFs."""
    out = tmp_path / "run_gwas"

    result = run_gwas(
        square_dataset["vcf"],
        square_dataset["phenotypes"],
        {"model": "linear", "trait": "yield", "quality_control": KEEP_ALL},
        output_dir=out,
    )

    assert result["status"] == "success"
    for row in read_table(out / "summary_statistics.tsv"):
        assert float(row["MAF"]) == pytest.approx(HAND_COUNTED_MAF[row["SNP"]], abs=5.1e-5), row["SNP"]
