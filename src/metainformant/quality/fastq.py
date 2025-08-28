"""FASTQ quality control analysis.

This module provides comprehensive FASTQ quality assessment similar to FastQC
but integrated into the METAINFORMANT framework. All implementations use real
computational methods without mocking.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union, Any
from collections import Counter, defaultdict
import warnings

import numpy as np
import pandas as pd

from ..core.io import open_text_auto
from ..dna.fastq import iter_fastq, FastqRecord


def analyze_fastq_quality(
    fastq_path: Union[str, Path],
    subsample: Optional[int] = None,
    seed: int = 42,
) -> Dict[str, Any]:
    """Comprehensive FASTQ quality analysis.
    
    Args:
        fastq_path: Path to FASTQ file (.fastq or .fastq.gz)
        subsample: Number of reads to analyze (None for all)
        seed: Random seed for subsampling
        
    Returns:
        Dictionary with comprehensive quality metrics
    """
    fastq_path = Path(fastq_path)
    if not fastq_path.exists():
        raise FileNotFoundError(f"FASTQ file not found: {fastq_path}")
    
    # Collect all reads for analysis
    reads = []
    total_reads = 0
    
    if subsample is not None:
        np.random.seed(seed)
    
    for read_id, sequence, quality in iter_fastq(fastq_path):
        total_reads += 1
        
        if subsample is None or len(reads) < subsample:
            reads.append(FastqRecord(read_id, sequence, quality))
        elif subsample is not None:
            # Reservoir sampling for large files
            replace_idx = np.random.randint(0, total_reads)
            if replace_idx < subsample:
                reads[replace_idx] = FastqRecord(read_id, sequence, quality)
    
    if len(reads) == 0:
        raise ValueError("No reads found in FASTQ file")
    
    print(f"Analyzing {len(reads)} reads from {fastq_path.name}")
    
    # Run all quality analyses
    results = {
        'basic_stats': basic_statistics(reads),
        'per_base_quality': per_base_quality(reads),
        'per_sequence_quality': per_sequence_quality(reads),
        'sequence_lengths': sequence_length_distribution(reads),
        'gc_content': gc_content_distribution(reads),
        'adapter_content': adapter_content(reads),
        'overrepresented_seqs': overrepresented_sequences(reads),
        'duplication_levels': duplication_levels(reads),
        'n_content': n_content_per_position(reads),
        'quality_scores': quality_score_distribution(reads),
    }
    
    # Add metadata
    results['metadata'] = {
        'filename': fastq_path.name,
        'total_reads_in_file': total_reads,
        'reads_analyzed': len(reads),
        'subsampled': subsample is not None,
    }
    
    return results


def basic_statistics(reads: List[FastqRecord]) -> Dict[str, Any]:
    """Calculate basic FASTQ statistics.
    
    Args:
        reads: List of FASTQ records
        
    Returns:
        Dictionary with basic statistics
    """
    if not reads:
        return {}
    
    lengths = [len(read.sequence) for read in reads]
    gc_counts = [read.sequence.count('G') + read.sequence.count('C') for read in reads]
    
    # Quality scores (convert from ASCII)
    quality_scores = []
    for read in reads:
        scores = [ord(char) - 33 for char in read.quality]  # Phred+33
        quality_scores.extend(scores)
    
    stats = {
        'total_sequences': len(reads),
        'sequence_length_min': min(lengths),
        'sequence_length_max': max(lengths),
        'sequence_length_mean': np.mean(lengths),
        'sequence_length_median': np.median(lengths),
        'gc_content_mean': np.mean([gc / length for gc, length in zip(gc_counts, lengths)]) * 100,
        'quality_score_mean': np.mean(quality_scores),
        'quality_score_median': np.median(quality_scores),
        'quality_score_min': min(quality_scores),
        'quality_score_max': max(quality_scores),
    }
    
    # Determine if sequences are uniform length
    stats['uniform_length'] = len(set(lengths)) == 1
    
    return stats


def per_base_quality(reads: List[FastqRecord]) -> Dict[str, Any]:
    """Calculate per-base quality statistics.
    
    Args:
        reads: List of FASTQ records
        
    Returns:
        Dictionary with per-base quality statistics
    """
    if not reads:
        return {}
    
    # Find maximum read length
    max_length = max(len(read.sequence) for read in reads)
    
    # Collect quality scores for each position
    position_qualities = defaultdict(list)
    
    for read in reads:
        for pos, qual_char in enumerate(read.quality):
            quality_score = ord(qual_char) - 33  # Phred+33
            position_qualities[pos].append(quality_score)
    
    # Calculate statistics per position
    per_base_stats = []
    for pos in range(max_length):
        if pos in position_qualities:
            scores = position_qualities[pos]
            stats = {
                'position': pos + 1,  # 1-based
                'mean': np.mean(scores),
                'median': np.median(scores),
                'q1': np.percentile(scores, 25),
                'q3': np.percentile(scores, 75),
                'min': min(scores),
                'max': max(scores),
                'count': len(scores),
            }
        else:
            # No data for this position
            stats = {
                'position': pos + 1,
                'mean': 0, 'median': 0, 'q1': 0, 'q3': 0,
                'min': 0, 'max': 0, 'count': 0,
            }
        per_base_stats.append(stats)
    
    return {
        'per_base_quality': per_base_stats,
        'max_position': max_length,
    }


def per_sequence_quality(reads: List[FastqRecord]) -> Dict[str, Any]:
    """Calculate per-sequence quality score distribution.
    
    Args:
        reads: List of FASTQ records
        
    Returns:
        Dictionary with per-sequence quality distribution
    """
    sequence_qualities = []
    
    for read in reads:
        # Calculate mean quality for this sequence
        quality_scores = [ord(char) - 33 for char in read.quality]
        mean_quality = np.mean(quality_scores)
        sequence_qualities.append(mean_quality)
    
    # Create histogram bins
    min_qual, max_qual = min(sequence_qualities), max(sequence_qualities)
    bins = np.linspace(min_qual, max_qual, 50)
    hist, bin_edges = np.histogram(sequence_qualities, bins=bins)
    
    return {
        'sequence_qualities': sequence_qualities,
        'histogram_counts': hist.tolist(),
        'histogram_bins': bin_edges.tolist(),
        'mean_quality_mean': np.mean(sequence_qualities),
        'mean_quality_median': np.median(sequence_qualities),
        'quality_distribution': dict(zip(
            [f"{edge:.1f}" for edge in bin_edges[:-1]], 
            hist.tolist()
        )),
    }


def sequence_length_distribution(reads: List[FastqRecord]) -> Dict[str, Any]:
    """Analyze sequence length distribution.
    
    Args:
        reads: List of FASTQ records
        
    Returns:
        Dictionary with length distribution data
    """
    lengths = [len(read.sequence) for read in reads]
    length_counts = Counter(lengths)
    
    return {
        'lengths': lengths,
        'length_distribution': dict(length_counts),
        'min_length': min(lengths),
        'max_length': max(lengths),
        'mean_length': np.mean(lengths),
        'median_length': np.median(lengths),
        'mode_length': max(length_counts, key=length_counts.get),
        'uniform_length': len(set(lengths)) == 1,
    }


def gc_content_distribution(reads: List[FastqRecord]) -> Dict[str, Any]:
    """Analyze GC content distribution.
    
    Args:
        reads: List of FASTQ records
        
    Returns:
        Dictionary with GC content analysis
    """
    gc_contents = []
    
    for read in reads:
        sequence = read.sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        at_count = sequence.count('A') + sequence.count('T')
        total_bases = gc_count + at_count
        
        if total_bases > 0:
            gc_percent = (gc_count / total_bases) * 100
        else:
            gc_percent = 0
        
        gc_contents.append(gc_percent)
    
    # Create histogram
    bins = np.arange(0, 101, 2)  # 0-100% in 2% bins
    hist, bin_edges = np.histogram(gc_contents, bins=bins)
    
    return {
        'gc_contents': gc_contents,
        'histogram_counts': hist.tolist(),
        'histogram_bins': bin_edges.tolist(),
        'mean_gc': np.mean(gc_contents),
        'median_gc': np.median(gc_contents),
        'gc_distribution': dict(zip(
            [f"{edge:.0f}" for edge in bin_edges[:-1]], 
            hist.tolist()
        )),
    }


def adapter_content(
    reads: List[FastqRecord],
    adapters: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Detect adapter contamination.
    
    Args:
        reads: List of FASTQ records
        adapters: Dictionary of adapter sequences {name: sequence}
        
    Returns:
        Dictionary with adapter contamination analysis
    """
    if adapters is None:
        # Common adapter sequences
        adapters = {
            'Illumina_Universal': 'AGATCGGAAGAG',
            'Illumina_Small_RNA': 'TGGAATTCTCGG',
            'Nextera_Transposase': 'CTGTCTCTTATA',
            'Poly_A': 'AAAAAAAAAA',
            'Poly_T': 'TTTTTTTTTT',
        }
    
    adapter_results = {}
    
    for adapter_name, adapter_seq in adapters.items():
        # Check for adapter presence at different positions
        position_contamination = defaultdict(int)
        total_contaminated = 0
        
        for read in reads:
            sequence = read.sequence.upper()
            adapter_upper = adapter_seq.upper()
            
            # Look for adapter starting at each position
            for start_pos in range(len(sequence) - len(adapter_seq) + 1):
                subseq = sequence[start_pos:start_pos + len(adapter_seq)]
                
                # Allow some mismatches (80% similarity)
                matches = sum(1 for a, b in zip(subseq, adapter_upper) if a == b)
                similarity = matches / len(adapter_seq)
                
                if similarity >= 0.8:  # 80% match threshold
                    position_contamination[start_pos] += 1
                    total_contaminated += 1
                    break  # Only count once per read
        
        adapter_results[adapter_name] = {
            'total_contaminated_reads': total_contaminated,
            'contamination_percentage': (total_contaminated / len(reads)) * 100,
            'position_contamination': dict(position_contamination),
        }
    
    return adapter_results


def overrepresented_sequences(
    reads: List[FastqRecord],
    min_count: int = 10,
    min_percentage: float = 0.1,
) -> Dict[str, Any]:
    """Identify overrepresented sequences.
    
    Args:
        reads: List of FASTQ records
        min_count: Minimum count to consider overrepresented
        min_percentage: Minimum percentage to consider overrepresented
        
    Returns:
        Dictionary with overrepresented sequence analysis
    """
    sequence_counts = Counter([read.sequence for read in reads])
    total_sequences = len(reads)
    
    overrepresented = {}
    for sequence, count in sequence_counts.most_common():
        percentage = (count / total_sequences) * 100
        
        if count >= min_count and percentage >= min_percentage:
            overrepresented[sequence] = {
                'count': count,
                'percentage': percentage,
                'length': len(sequence),
                'possible_source': _identify_sequence_source(sequence),
            }
        
        # Limit to top 20 overrepresented sequences
        if len(overrepresented) >= 20:
            break
    
    return {
        'overrepresented_sequences': overrepresented,
        'total_unique_sequences': len(sequence_counts),
        'most_common_count': sequence_counts.most_common(1)[0][1] if sequence_counts else 0,
    }


def duplication_levels(reads: List[FastqRecord]) -> Dict[str, Any]:
    """Analyze sequence duplication levels.
    
    Args:
        reads: List of FASTQ records
        
    Returns:
        Dictionary with duplication analysis
    """
    sequence_counts = Counter([read.sequence for read in reads])
    
    # Count how many sequences appear at each duplication level
    duplication_counts = Counter(sequence_counts.values())
    
    # Calculate statistics
    total_sequences = len(reads)
    unique_sequences = len(sequence_counts)
    
    duplication_levels_data = {}
    for dup_level, count in duplication_counts.items():
        duplication_levels_data[dup_level] = {
            'sequences_with_this_dup_level': count,
            'percentage_of_unique': (count / unique_sequences) * 100,
            'total_reads_from_this_level': dup_level * count,
            'percentage_of_total': (dup_level * count / total_sequences) * 100,
        }
    
    return {
        'duplication_levels': duplication_levels_data,
        'total_sequences': total_sequences,
        'unique_sequences': unique_sequences,
        'duplication_rate': (total_sequences - unique_sequences) / total_sequences * 100,
        'mean_duplication': np.mean(list(sequence_counts.values())),
    }


def n_content_per_position(reads: List[FastqRecord]) -> Dict[str, Any]:
    """Analyze N content at each position.
    
    Args:
        reads: List of FASTQ records
        
    Returns:
        Dictionary with N content per position
    """
    if not reads:
        return {}
    
    max_length = max(len(read.sequence) for read in reads)
    position_n_counts = defaultdict(int)
    position_totals = defaultdict(int)
    
    for read in reads:
        for pos, base in enumerate(read.sequence.upper()):
            position_totals[pos] += 1
            if base == 'N':
                position_n_counts[pos] += 1
    
    per_position_n = []
    for pos in range(max_length):
        n_count = position_n_counts.get(pos, 0)
        total = position_totals.get(pos, 0)
        percentage = (n_count / total * 100) if total > 0 else 0
        
        per_position_n.append({
            'position': pos + 1,
            'n_count': n_count,
            'total_bases': total,
            'n_percentage': percentage,
        })
    
    return {
        'per_position_n_content': per_position_n,
        'overall_n_percentage': sum(position_n_counts.values()) / sum(position_totals.values()) * 100,
    }


def quality_score_distribution(reads: List[FastqRecord]) -> Dict[str, Any]:
    """Analyze quality score distribution.
    
    Args:
        reads: List of FASTQ records
        
    Returns:
        Dictionary with quality score distribution
    """
    all_quality_scores = []
    
    for read in reads:
        quality_scores = [ord(char) - 33 for char in read.quality]  # Phred+33
        all_quality_scores.extend(quality_scores)
    
    if not all_quality_scores:
        return {}
    
    # Create histogram
    bins = np.arange(0, 42)  # Quality scores typically 0-41
    hist, bin_edges = np.histogram(all_quality_scores, bins=bins)
    
    return {
        'quality_scores': all_quality_scores,
        'histogram_counts': hist.tolist(),
        'histogram_bins': bin_edges.tolist(),
        'mean_quality': np.mean(all_quality_scores),
        'median_quality': np.median(all_quality_scores),
        'quality_distribution': dict(zip(
            [int(edge) for edge in bin_edges[:-1]], 
            hist.tolist()
        )),
    }


def _identify_sequence_source(sequence: str) -> str:
    """Identify possible source of overrepresented sequence.
    
    Args:
        sequence: DNA sequence
        
    Returns:
        Possible source description
    """
    sequence = sequence.upper()
    
    # Check for common patterns
    if sequence == 'A' * len(sequence):
        return 'Poly-A tail'
    elif sequence == 'T' * len(sequence):
        return 'Poly-T sequence'
    elif sequence == 'G' * len(sequence):
        return 'Poly-G sequence'
    elif sequence == 'C' * len(sequence):
        return 'Poly-C sequence'
    elif 'AGATCGGAAGAG' in sequence:
        return 'Illumina adapter'
    elif 'TGGAATTCTCGG' in sequence:
        return 'Illumina small RNA adapter'
    elif 'CTGTCTCTTATA' in sequence:
        return 'Nextera adapter'
    elif len(set(sequence)) <= 2:
        return 'Low complexity sequence'
    elif sequence.count('CG') / len(sequence) > 0.3:
        return 'High CG content (possible ribosomal)'
    else:
        return 'Unknown'
