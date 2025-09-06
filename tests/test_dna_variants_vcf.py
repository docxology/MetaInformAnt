from __future__ import annotations

from pathlib import Path

from metainformant.dna import variants


def test_parse_vcf_minimal(tmp_path: Path) -> None:
    vcf_text = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t1\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1
chr1\t2\trs2\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/1
"""
    p = tmp_path / "small.vcf"
    p.write_text(vcf_text)
    parsed = variants.parse_vcf(p)
    assert parsed["num_variants"] == 2
    assert parsed["samples"] == ["S1", "S2"]
