"""Shim module for metainformant.dna.io.fasta.

Provides backward compatible imports for FASTQ handling.
"""

from .fastq import FastqRecord, read_fastq, write_fastq, assess_quality, filter_reads, convert_fastq_to_fasta, trim_reads, summarize_fastq, calculate_per_base_quality, iter_fastq, gc_content, average_phred_by_position
