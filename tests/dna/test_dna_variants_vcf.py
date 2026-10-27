"""Tests for DNA variant VCF parsing."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from metainformant.dna.variation import variants


def test_parse_vcf_minimal(tmp_path: Path) -> None:
    """Test parsing a minimal VCF file."""
    vcf_text = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t1\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1
chr1\t2\trs2\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/1
"""
    p = tmp_path / "small.vcf"
    p.write_text(vcf_text)
    parsed = variants.parse_vcf(p)
    # Function returns 'total_variants' not 'num_variants'
    assert parsed["total_variants"] == 2
    assert parsed["samples"] == ["S1", "S2"]
    assert "variants" in parsed


def test_parse_vcf_returns_dict(tmp_path: Path) -> None:
    """Test that parse_vcf returns a dictionary."""
    vcf_text = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1\t.\tA\tG\t.\tPASS\t.
"""
    p = tmp_path / "minimal.vcf"
    p.write_text(vcf_text)
    parsed = variants.parse_vcf(p)
    assert isinstance(parsed, dict)


VCF_TEXT = (
    "##fileformat=VCFv4.2\n"
    "##source=metainformant-test\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
    "chr1\t100\trs1\tA\tG\t50.5\tPASS\tDP=14;AF=0.5;DB\tGT:DP\t0/1:7\t1/1:7\n"
    "chr1\t250\t.\tT\tC,G\t.\tq10\tDP=9\tGT:DP\t0/0:4\t0/2:5\n"
    "chr2\t42\trs3\tGA\tG\t99\tPASS\t.\tGT:DP\t1/1:12\t./.:0\n"
)


def test_parse_vcf_reads_gzip_compressed_files(tmp_path: Path) -> None:
    """A ``.vcf.gz`` used to raise ``UnicodeDecodeError`` (parse_vcf_full has always read it)."""
    plain = tmp_path / "calls.vcf"
    plain.write_text(VCF_TEXT)
    compressed = tmp_path / "calls.vcf.gz"
    with gzip.open(compressed, "wt") as handle:
        handle.write(VCF_TEXT)

    parsed = variants.parse_vcf(compressed)

    assert parsed["total_variants"] == 3
    assert parsed["samples"] == ["S1", "S2"]
    assert parsed["metadata"]["fileformat"] == "VCFv4.2"
    first, second, third = parsed["variants"]
    assert (first["chrom"], first["pos"], first["id"], first["ref"], first["alt"]) == ("chr1", 100, "rs1", "A", ["G"])
    assert first["qual"] == 50.5
    assert first["info"] == {"DP": 14, "AF": 0.5, "DB": True}
    assert first["samples"]["S1"] == {"GT": "0/1", "DP": "7"}
    assert second["alt"] == ["C", "G"] and second["qual"] is None and second["filter"] == ["q10"]
    assert third["ref"] == "GA" and third["samples"]["S2"] == {"GT": "./.", "DP": "0"}
    # identical to what the plain-text file gives
    assert parsed == variants.parse_vcf(plain)


def test_parse_vcf_reads_concatenated_gzip_members(tmp_path: Path) -> None:
    """bgzip writes a series of gzip members; every one of them must be read, not just the first."""
    lines = VCF_TEXT.splitlines(keepends=True)
    members = [gzip.compress("".join(lines[:4]).encode()), gzip.compress("".join(lines[4:]).encode())]
    compressed = tmp_path / "blocks.vcf.gz"
    compressed.write_bytes(b"".join(members))
    plain = tmp_path / "blocks.vcf"
    plain.write_text(VCF_TEXT)

    parsed = variants.parse_vcf(compressed)

    assert parsed["total_variants"] == 3
    assert parsed == variants.parse_vcf(plain)


def test_parse_vcf_missing_gzip_file_still_raises_file_not_found(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        variants.parse_vcf(tmp_path / "missing.vcf.gz")


def test_write_vcf_to_a_gz_path_is_real_gzip_that_parse_vcf_reads_back(tmp_path: Path) -> None:
    """The writer must agree with the reader about what a ``.gz`` name means."""
    plain = tmp_path / "in.vcf"
    plain.write_text(VCF_TEXT)
    data = variants.parse_vcf(plain)

    written_plain = tmp_path / "out.vcf"
    written_gz = tmp_path / "out.vcf.gz"
    variants.write_vcf(data, written_plain)
    variants.write_vcf(data, written_gz)

    assert written_gz.read_bytes()[:2] == b"\x1f\x8b"  # gzip magic, not plain text under a .gz name
    assert gzip.decompress(written_gz.read_bytes()).decode() == written_plain.read_text()
    assert variants.parse_vcf(written_gz) == variants.parse_vcf(written_plain)
    assert variants.parse_vcf(written_gz)["total_variants"] == 3
