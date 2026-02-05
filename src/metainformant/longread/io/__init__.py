"""Long-read I/O module for FAST5, POD5, BAM, and format conversion.

Handles reading and writing of long-read sequencing data formats including
ONT FAST5/POD5 signal data, long-read BAM files with methylation tags,
and format conversions between common long-read formats.
"""

from __future__ import annotations

from . import fast5
from . import bam
from . import formats

__all__ = [
    "fast5",
    "bam",
    "formats",
]
