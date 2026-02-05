"""Shim module for metainformant.dna.io.fasta.

Provides backward compatible imports for FASTQ handling.
"""

from .fastq import (
    FastqRecord,
    assess_quality,
    average_phred_by_position,
    calculate_per_base_quality,
    convert_fastq_to_fasta,
    filter_reads,
    gc_content,
    iter_fastq,
    read_fastq,
    summarize_fastq,
    trim_reads,
    write_fastq,
)
