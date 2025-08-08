from __future__ import annotations

from pathlib import Path
from typing import Any, Dict


def parse_vcf(path: str | Path) -> Dict[str, Any]:
    """Parse a minimal VCF to extract sample names and variant count.

    Only inspects header for samples and counts non-header lines as variants.
    """
    samples: list[str] = []
    num_variants = 0
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                samples = parts[9:]
                continue
            if not line.strip():
                continue
            num_variants += 1
    return {"samples": samples, "num_variants": num_variants}


