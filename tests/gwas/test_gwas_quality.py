"""Tests for GWAS quality control functions."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest

from metainformant.gwas.analysis.quality import apply_qc_filters, parse_vcf_full, write_filtered_vcf


def test_parse_vcf_full_basic(tmp_path: Path) -> None:
    """Test parsing VCF file with genotypes."""
    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	rs1	A	G	60	PASS	DP=30	GT	0/1	1/1	0/0
chr1	200	rs2	T	C	80	PASS	DP=40	GT	0/0	0/1	1/1
chr1	300	rs3	G	A	50	PASS	DP=20	GT	1/1	0/0	0/1
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    result = parse_vcf_full(vcf_file)

    assert result["samples"] == ["S1", "S2", "S3"]
    assert len(result["variants"]) == 3
    assert len(result["genotypes"]) == 3  # 3 samples
    assert len(result["genotypes"][0]) == 3  # 3 variants

    # Check first variant genotypes: S1=0/1(1), S2=1/1(2), S3=0/0(0)
    assert result["genotypes"][0][0] == 1  # S1, variant 0
    assert result["genotypes"][1][0] == 2  # S2, variant 0
    assert result["genotypes"][2][0] == 0  # S3, variant 0


def test_parse_vcf_full_uses_safe_sidecar_cache(tmp_path: Path) -> None:
    """Second parse should use the pickle-free cache when source metadata matches."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	rs1	A	G	60	PASS	.	GT	0/1
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    first = parse_vcf_full(vcf_file)
    cache_file = tmp_path / "test.vcf.npz"
    assert cache_file.exists()

    second = parse_vcf_full(vcf_file)
    assert second["samples"] == first["samples"]
    assert second["variants"] == first["variants"]
    assert second["metadata"] == first["metadata"]
    assert second["genotypes"].tolist() == first["genotypes"].tolist()


def test_parse_vcf_full_missing_genotypes(tmp_path: Path) -> None:
    """Test parsing VCF with missing genotypes."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	rs1	A	G	60	PASS	.	GT	0/1	./.	0/0
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    result = parse_vcf_full(vcf_file)
    assert len(result["genotypes"]) == 3
    assert result["genotypes"][0][0] == 1  # S1: 0/1 = 1
    assert result["genotypes"][1][0] == -1  # S2: ./. = missing
    assert result["genotypes"][2][0] == 0  # S3: 0/0 = 0


def test_apply_qc_filters_basic(tmp_path: Path) -> None:
    """Test applying basic QC filters."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
chr1	100	rs1	A	G	60	PASS	.	GT	0/1	1/1	0/0	0/1	0/0
chr1	200	rs2	T	C	80	PASS	.	GT	0/0	0/0	0/0	0/0	0/0
chr1	300	rs3	G	A	50	PASS	.	GT	1/1	1/1	1/1	1/1	1/1
chr1	400	rs4	C	T	30	PASS	.	GT	./.	./.	./.	./.	./.
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    qc_config = {
        "min_maf": 0.1,
        "max_missing": 0.2,
        "min_qual": 40.0,
        "exclude_indels": False,
    }

    result = apply_qc_filters(vcf_file, qc_config)

    assert result["status"] == "success"
    assert result["num_variants_before"] == 4
    assert result["num_variants_after"] <= 4

    # rs2 should be filtered (all hom ref = MAF = 0)
    # rs3 should be filtered (all hom alt = MAF = 0)
    # rs4 should be filtered (all missing)
    # rs1 should pass (MAF = 0.4, qual = 60)


def test_apply_qc_filters_maf(tmp_path: Path) -> None:
    """Test MAF filtering."""
    # Create VCF with variants at different MAF levels
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10
chr1	100	rs1	A	G	60	PASS	.	GT	0/1	1/1	0/0	0/1	0/0	0/1	0/0	0/1	0/0	0/1
chr1	200	rs2	T	C	80	PASS	.	GT	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	0/0	1/1
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    # rs1: 5 het + 1 hom alt = 6 alt alleles out of 20 = MAF = 0.3
    # rs2: 1 hom alt = 2 alt alleles out of 20 = MAF = 0.1

    qc_config = {"min_maf": 0.15, "max_missing": 1.0, "min_qual": 0.0}
    result = apply_qc_filters(vcf_file, qc_config)

    assert result["status"] == "success"
    # Only rs1 should pass (MAF >= 0.15)
    assert 1 <= result["num_variants_after"] <= 2


def test_apply_qc_filters_missing_data(tmp_path: Path) -> None:
    """Test missing data filtering."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
chr1	100	rs1	A	G	60	PASS	.	GT	0/1	1/1	0/0	0/1	0/0
chr1	200	rs2	T	C	80	PASS	.	GT	0/0	./.	./.	./.	./.
chr1	300	rs3	G	A	50	PASS	.	GT	1/1	1/1	./.	./.	./.
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    qc_config = {"min_maf": 0.0, "max_missing": 0.5, "min_qual": 0.0}
    result = apply_qc_filters(vcf_file, qc_config)

    assert result["status"] == "success"
    # rs2 has 80% missing (4/5), rs3 has 60% missing (3/5)
    # Only rs1 should pass (0% missing)


def test_apply_qc_filters_hwe(tmp_path: Path) -> None:
    """Test Hardy-Weinberg equilibrium filtering."""
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10
chr1	100	rs1	A	G	60	PASS	.	GT	0/0	0/0	0/0	0/0	0/0	1/1	1/1	1/1	1/1	1/1
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    # This variant violates HWE (no heterozygotes, extreme deviation)
    qc_config = {"min_maf": 0.0, "max_missing": 1.0, "min_qual": 0.0, "hwe_pval": 1e-6}
    result = apply_qc_filters(vcf_file, qc_config)

    assert result["status"] == "success"
    # Variant with extreme HWE violation might be filtered


QC_ROUNDTRIP_VCF = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t100\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t1/1\t0/0\t0/1
chr1\t200\trs2\tA\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1\t0/0
chr1\t300\trs3\tA\tC\t50\tPASS\tDP=10\tGT\t0/0\t0/1\t1/1
"""


def test_write_filtered_vcf_round_trip_preserves_orientation(tmp_path: Path) -> None:
    """Written genotypes must keep sample-major orientation: column i per variant."""
    vcf_file = tmp_path / "in.vcf"
    vcf_file.write_text(QC_ROUNDTRIP_VCF)
    parsed = parse_vcf_full(vcf_file)

    out_vcf = tmp_path / "out.vcf"
    write_filtered_vcf(parsed, out_vcf)

    lines = out_vcf.read_text().splitlines()
    data_lines = [ln for ln in lines if not ln.startswith("#")]
    # rs1 columns are S1=1/1, S2=0/0, S3=0/1 — sample order must be preserved
    assert data_lines[0].split("\t")[9:] == ["1/1", "0/0", "0/1"]
    assert data_lines[1].split("\t")[9:] == ["0/1", "1/1", "0/0"]
    assert data_lines[2].split("\t")[9:] == ["0/0", "0/1", "1/1"]

    reparsed = parse_vcf_full(out_vcf)
    assert reparsed["samples"] == ["S1", "S2", "S3"]
    assert reparsed["genotypes"].tolist() == [[2, 1, 0], [0, 2, 1], [1, 0, 2]]


def test_write_filtered_vcf_after_qc_keeps_sample_columns(tmp_path: Path) -> None:
    """QC output written back to VCF must retain per-sample genotype columns."""
    vcf_file = tmp_path / "in.vcf"
    vcf_file.write_text(QC_ROUNDTRIP_VCF)
    qc = apply_qc_filters(vcf_file, {"min_maf": 0.0, "max_missing": 1.0, "min_hwe_p": 0.0})
    assert qc["num_variants_after"] == 3

    out_vcf = tmp_path / "qc.vcf"
    write_filtered_vcf(qc["filtered_data"], out_vcf)

    reparsed = parse_vcf_full(out_vcf)
    assert reparsed["genotypes"].tolist() == [[2, 1, 0], [0, 2, 1], [1, 0, 2]]


# ---------------------------------------------------------------------------
# QUAL handling: min_qual filter and the QUAL column of written VCFs
# ---------------------------------------------------------------------------

QUAL_VCF = (
    "##fileformat=VCFv4.2\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
    "chr1\t100\trsLOW\tA\tG\t10\tPASS\t.\tGT\t0/1\t0/0\t0/1\t1/1\n"
    "chr1\t200\trsHIGH\tT\tC\t60.5\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\n"
    "chr1\t300\trsMISS\tG\tA\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\t0/1\n"
    "chr1\t400\trsEDGE\tC\tT\t30\tPASS\t.\tGT\t1/1\t0/1\t0/0\t0/1\n"
    "chr1\t500\trsTOP\tA\tC\t999\tPASS\t.\tGT\t0/1\t0/0\t1/1\t0/1\n"
)

# Every variant is polymorphic and complete, so only QUAL can remove a variant under this configuration.
QC_PASS_THROUGH = {"min_maf": 0.0, "max_missing": 1.0, "min_hwe_p": 0.0}

_GT_TO_DOSAGE = {"0/0": 0, "0/1": 1, "1/1": 2, "./.": -1}


def _raw_vcf_records(path: Path) -> list[dict[str, Any]]:
    """Read VCF data lines with plain string handling, independent of parse_vcf_full."""
    records: list[dict[str, Any]] = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        records.append({"id": columns[2], "qual": columns[5], "gt": columns[9:]})
    return records


def test_apply_qc_filters_min_qual_removes_low_and_missing_qual(tmp_path: Path) -> None:
    """min_qual keeps QUAL >= threshold (inclusive) and removes low and missing QUAL."""
    vcf_file = tmp_path / "qual.vcf"
    vcf_file.write_text(QUAL_VCF)
    records = _raw_vcf_records(vcf_file)

    # Oracle: plain string parsing of the QUAL column. A missing QUAL ('.') cannot meet a threshold.
    expected_ids = [r["id"] for r in records if r["qual"] != "." and float(r["qual"]) >= 30.0]
    assert expected_ids == ["rsHIGH", "rsEDGE", "rsTOP"]

    unfiltered = apply_qc_filters(vcf_file, dict(QC_PASS_THROUGH))
    assert unfiltered["num_variants_after"] == 5

    filtered = apply_qc_filters(vcf_file, {**QC_PASS_THROUGH, "min_qual": 30.0})

    assert [v["id"] for v in filtered["filtered_data"]["variants"]] == expected_ids
    assert filtered["num_variants_before"] == 5
    assert filtered["num_variants_after"] == 3
    assert filtered["qc_filters"]["min_qual"] == 30.0
    assert filtered["qual_filter"] == {"min_qual": 30.0, "n_removed": 2, "n_missing_qual": 1}

    # Genotype columns must stay aligned with the surviving variants (sample-major layout).
    kept = [r for r in records if r["id"] in expected_ids]
    expected_genotypes = [[_GT_TO_DOSAGE[r["gt"][s]] for r in kept] for s in range(4)]
    assert np.asarray(filtered["filtered_data"]["genotypes"]).tolist() == expected_genotypes

    # A pre-parsed dictionary takes the same path as a file path.
    from_parsed = apply_qc_filters(parse_vcf_full(vcf_file), {**QC_PASS_THROUGH, "min_qual": 30.0})
    assert [v["id"] for v in from_parsed["filtered_data"]["variants"]] == expected_ids


def test_apply_qc_filters_min_qual_composes_with_other_filters(tmp_path: Path) -> None:
    """Indices from the QUAL stage must map back correctly through the MAF stage."""
    vcf_text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
        "chr1\t100\trsLOWQ\tA\tG\t5\tPASS\t.\tGT\t0/1\t0/0\t0/1\t1/1\n"
        "chr1\t200\trsMONO\tT\tC\t80\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\n"
        "chr1\t300\trsGOOD\tG\tA\t80\tPASS\t.\tGT\t0/1\t0/0\t1/1\t0/1\n"
    )
    vcf_file = tmp_path / "compose.vcf"
    vcf_file.write_text(vcf_text)

    result = apply_qc_filters(vcf_file, {"min_maf": 0.1, "max_missing": 1.0, "min_hwe_p": 0.0, "min_qual": 30.0})

    # rsLOWQ fails QUAL, rsMONO fails MAF; only rsGOOD (the third input variant) survives.
    assert [v["id"] for v in result["filtered_data"]["variants"]] == ["rsGOOD"]
    assert np.asarray(result["filtered_data"]["genotypes"]).tolist() == [[1], [0], [2], [1]]


def _with_all_qual_missing(vcf_text: str) -> str:
    """Return the VCF with every data line's QUAL column replaced by '.'."""
    out: list[str] = []
    for line in vcf_text.splitlines(keepends=True):
        if line.startswith("#"):
            out.append(line)
            continue
        columns = line.rstrip("\n").split("\t")
        columns[5] = "."
        out.append("\t".join(columns) + "\n")
    return "".join(out)


def test_apply_qc_filters_min_qual_with_all_qual_missing_removes_everything(tmp_path: Path) -> None:
    """The missing-QUAL policy is uniform: with a threshold set, '.' never passes (even for 0)."""
    vcf_file = tmp_path / "noqual.vcf"
    vcf_file.write_text(_with_all_qual_missing(QUAL_VCF))
    assert {r["qual"] for r in _raw_vcf_records(vcf_file)} == {"."}

    without_threshold = apply_qc_filters(vcf_file, dict(QC_PASS_THROUGH))
    assert without_threshold["num_variants_after"] == 5

    with_threshold = apply_qc_filters(vcf_file, {**QC_PASS_THROUGH, "min_qual": 0.0})
    assert with_threshold["status"] == "success"
    assert with_threshold["num_variants_after"] == 0
    assert with_threshold["qual_filter"]["n_missing_qual"] == 5
    assert with_threshold["qual_filter"]["n_removed"] == 5


def test_apply_qc_filters_min_qual_input_validation(tmp_path: Path) -> None:
    """A malformed threshold or QUAL-less variant data is an error, not a silent no-op."""
    vcf_file = tmp_path / "qual.vcf"
    vcf_file.write_text(QUAL_VCF)

    with pytest.raises(ValueError, match="min_qual"):
        apply_qc_filters(vcf_file, {**QC_PASS_THROUGH, "min_qual": "high"})
    with pytest.raises(ValueError, match="min_qual"):
        apply_qc_filters(vcf_file, {**QC_PASS_THROUGH, "min_qual": float("nan")})

    parsed = parse_vcf_full(vcf_file)
    parsed["variants"] = parsed["variants"][:2]  # fewer variant records than genotype columns
    with pytest.raises(ValueError, match="variant"):
        apply_qc_filters(parsed, {**QC_PASS_THROUGH, "min_qual": 30.0})


def test_write_filtered_vcf_writes_dot_for_missing_qual(tmp_path: Path) -> None:
    """A missing QUAL is written as '.', never as the text 'None'."""
    vcf_file = tmp_path / "qual.vcf"
    vcf_file.write_text(QUAL_VCF)
    parsed = parse_vcf_full(vcf_file)
    assert [v["qual"] for v in parsed["variants"]] == [10.0, 60.5, None, 30.0, 999.0]

    out_vcf = tmp_path / "out.vcf"
    write_filtered_vcf(parsed, out_vcf)

    assert "None" not in out_vcf.read_text()
    written = [r["qual"] for r in _raw_vcf_records(out_vcf)]
    assert written[2] == "."
    assert [float(q) for q in written if q != "."] == [10.0, 60.5, 30.0, 999.0]

    reparsed = parse_vcf_full(out_vcf)
    assert [v["qual"] for v in reparsed["variants"]] == [10.0, 60.5, None, 30.0, 999.0]


def test_write_filtered_vcf_renders_none_fields_as_vcf_missing_marker(tmp_path: Path) -> None:
    """None/NaN in QUAL, ID, FILTER, INFO and FORMAT become VCF-valid text, not 'None'/'nan'."""
    data = {
        "samples": ["S1", "S2"],
        "metadata": {},
        "genotypes": [[1], [0]],  # sample-major: S1 and S2, one variant
        "variants": [
            {
                "chrom": "chr1",
                "pos": 10,
                "id": None,
                "ref": "A",
                "alt": ["G"],
                "qual": float("nan"),
                "filter": None,
                "info": None,
                "format": None,
            }
        ],
    }
    out_vcf = tmp_path / "none_fields.vcf"
    write_filtered_vcf(data, out_vcf)

    data_lines = [ln for ln in out_vcf.read_text().splitlines() if not ln.startswith("#")]
    assert data_lines[0].split("\t") == ["chr1", "10", ".", "A", "G", ".", ".", ".", "GT", "0/1", "0/0"]
