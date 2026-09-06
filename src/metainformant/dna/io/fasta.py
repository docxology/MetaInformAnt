"""Compatibility shim module for ``metainformant.dna.io.fasta``.

This module historically re-exported FASTA helpers. FASTA reading/writing now
lives in :mod:`metainformant.dna.sequence.core` (``read_fasta`` /
``write_fasta``); the module is kept so existing
``metainformant.dna.io.fasta`` imports continue to resolve.
"""
