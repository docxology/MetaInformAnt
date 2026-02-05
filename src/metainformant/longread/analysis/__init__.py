"""Long-read analysis module for modified bases, structural variants, and phasing.

Includes modified base detection (5mC, 6mA), structural variant calling from
split/supplementary alignments, and haplotype phasing using heterozygous variants.
"""

from __future__ import annotations

from . import modified_bases
from . import structural
from . import phasing

__all__ = [
    "modified_bases",
    "structural",
    "phasing",
]
