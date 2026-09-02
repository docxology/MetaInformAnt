"""Shared zero-mocks test support: real tiny VCF built with reportlab-free tools.

Builds a minimal but valid compressed+indexed VCF using real pysam/bcftools
paths where available. When bcftools is present, we construct an uncompressed
VCF text and compress/index it with actual bcftools calls.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

TINY_VCF_TEXT = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t1000\trs1\tA\tG\t50\tPASS\tAF=0.33\tGT\t0/0\t0/1\t0/1
chr1\t2000\trs2\tT\tC\t60\tPASS\tAF=0.5\tGT\t0/1\t0/1\t1/1
chr1\t3000\trs3\tG\tA\t45\tPASS\tAF=0.17\tGT\t0/0\t0/0\t0/1
"""


def build_tiny_vcf(tmp_path: Path) -> Path:
    """Write, compress, and index a tiny real VCF using bcftools.

    Args:
        tmp_path: pytest tmp_path; a ``vcf`` subdir is created.

    Returns:
        Path to the compressed VCF (``tiny.vcf.gz``).
    """
    vcf_dir = tmp_path / "vcf"
    vcf_dir.mkdir(parents=True, exist_ok=True)
    plain = vcf_dir / "tiny.vcf"
    plain.write_text(TINY_VCF_TEXT)
    compressed = vcf_dir / "tiny.vcf.gz"
    subprocess.run(
        ["bcftools", "view", "-Oz", "-o", str(compressed), str(plain)],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        ["bcftools", "index", str(compressed)], check=True, capture_output=True
    )
    return compressed
