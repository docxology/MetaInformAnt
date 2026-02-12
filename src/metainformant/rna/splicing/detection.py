"""Alternative splicing detection and quantification for RNA-seq data.

This module provides tools for identifying splice junctions from aligned
reads, classifying alternative splicing events, computing Percent Spliced
In (PSI) values, testing for differential splicing between conditions,
discovering novel junctions, and scoring splice site strength.

All implementations are pure Python with optional numpy/scipy acceleration.
No external R or specialized splicing tool dependencies required.

Main Functions:
    Junction Detection:
        - detect_splice_junctions: Identify splice junctions from aligned reads
        - find_novel_junctions: Discover unannotated splice junctions

    Event Classification:
        - classify_splicing_events: Classify junctions by splicing event type

    Quantification:
        - compute_psi: Percent Spliced In with confidence interval
        - differential_splicing: Test for differential splicing between groups

    Sequence Analysis:
        - compute_splice_site_strength: Score splice site using position weight matrix

Example:
    >>> from metainformant.rna.splicing import detection
    >>> alignments = [
    ...     {"chrom": "chr1", "start": 100, "end": 200, "cigar": "50M100N50M"},
    ...     {"chrom": "chr1", "start": 100, "end": 200, "cigar": "50M100N50M"},
    ... ]
    >>> junctions = detection.detect_splice_junctions(alignments, min_reads=1)
    >>> psi = detection.compute_psi(inclusion_reads=30, exclusion_reads=10)

.. deprecated::
    This module is a backward-compatibility shim.  The implementation has
    been split into :mod:`splice_sites` and :mod:`splice_analysis`.
    Import directly from those modules for new code.
"""

from __future__ import annotations

# Re-export everything from the split modules so that existing imports
# like ``from metainformant.rna.splicing.detection import detect_splice_junctions``
# continue to work.

from .splice_sites import (
    JunctionType,
    _ACCEPTOR_PWM,
    _DONOR_PWM,
    _classify_junction_type,
    _parse_cigar_junctions,
    compute_splice_site_strength,
    detect_splice_junctions,
)
from .splice_analysis import (
    SplicingEventType,
    _beta_confidence_interval,
    _classify_by_coordinates,
    _classify_with_gene_model,
    _empirical_bayes_test,
    _norm_ppf,
    _normal_cdf,
    _permutation_test,
    _wilcoxon_test,
    classify_splicing_events,
    compute_psi,
    differential_splicing,
    find_novel_junctions,
)

__all__ = [
    # Types
    "SplicingEventType",
    "JunctionType",
    # Junction detection (splice_sites)
    "detect_splice_junctions",
    "compute_splice_site_strength",
    # Event classification & quantification (splice_analysis)
    "classify_splicing_events",
    "compute_psi",
    "differential_splicing",
    "find_novel_junctions",
]
