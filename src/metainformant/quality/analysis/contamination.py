"""Contamination detection and analysis.

This module provides algorithms for detecting various types of contamination
in biological sequencing data, including microbial contamination, cross-species
contamination, and technical artifacts.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

from metainformant.core.data import validation
from metainformant.core.utils import errors
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


class ContaminationDetector:
    """Detector for various types of contamination in sequencing data."""

    def __init__(self, reference_genomes: Optional[Dict[str, str]] = None):
        """Initialize the contamination detector.

        Args:
            reference_genomes: Dictionary mapping species names to reference genome sequences
        """
        self.reference_genomes = reference_genomes or {}
        self.known_contaminants = self._load_known_contaminants()
        logger.info(f"Initialized contamination detector with {len(self.reference_genomes)} reference genomes")

    def _load_known_contaminants(self) -> Dict[str, str]:
        """Load known contaminant sequences."""
        # Common laboratory contaminants
        contaminants = {
            "ecoli": "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC",
            "yeast": "TCTGGACCTGCGGCTTAATTTGACTCAACACGGGAAACCTCACCCGGTTTGCTGGGTCGAGTTGCAGCCTTTCATCGCTGTGAAGCAGACCTTCGGCGGTGCGGTGGTGTAGGCCTGGGGGTTGGGCGGGGGGCTCGGGGCGGGG",
            "human": "CCTGCAGGTTGAAGCGGGAAGAGTGGAGGTTGCCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCG",
        }
        return contaminants

    def detect_microbial_contamination(self, sequences: List[str], threshold: float = 0.01) -> Dict[str, Any]:
        """Detect microbial contamination in sequences.

        Args:
            sequences: List of DNA sequences to analyze
            threshold: Minimum fraction of reads that must match a contaminant

        Returns:
            Dictionary with contamination detection results
        """
        if not sequences:
            return {"detected": False, "contaminants": []}

        logger.info(f"Detecting microbial contamination in {len(sequences)} sequences")

        results = []
        total_sequences = len(sequences)

        for contaminant_name, contaminant_seq in self.known_contaminants.items():
            matches = 0
            match_positions = []

            # Use k-mer based detection for efficiency
            k = 20  # k-mer size
            contaminant_kmers = set()
            for i in range(len(contaminant_seq) - k + 1):
                contaminant_kmers.add(contaminant_seq[i : i + k])

            for seq in sequences:
                seq_kmers = set()
                for i in range(len(seq) - k + 1):
                    seq_kmers.add(seq[i : i + k])

                # Calculate Jaccard similarity
                intersection = len(contaminant_kmers.intersection(seq_kmers))
                union = len(contaminant_kmers.union(seq_kmers))

                if union > 0:
                    similarity = intersection / union
                    if similarity > 0.1:  # Significant k-mer overlap
                        matches += 1
                        match_positions.append(similarity)

            contamination_rate = matches / total_sequences

            if contamination_rate >= threshold:
                results.append(
                    {
                        "contaminant": contaminant_name,
                        "contamination_rate": contamination_rate,
                        "matches": matches,
                        "average_similarity": np.mean(match_positions) if match_positions else 0,
                    }
                )

        detected = len(results) > 0

        return {
            "detected": detected,
            "contaminants": results,
            "total_sequences_analyzed": total_sequences,
            "threshold": threshold,
        }

    def detect_cross_species_contamination(
        self, sequences: List[str], target_species: str, other_species: List[str]
    ) -> Dict[str, Any]:
        """Detect cross-species contamination.

        Args:
            sequences: List of sequences to analyze
            target_species: Expected species for the sequences
            other_species: List of potential contaminant species

        Returns:
            Dictionary with cross-species contamination results
        """
        if not sequences or target_species not in self.reference_genomes:
            return {"detected": False, "contaminants": []}

        logger.info(f"Detecting cross-species contamination (target: {target_species})")

        target_genome = self.reference_genomes[target_species]
        results = []

        for species in other_species:
            if species not in self.reference_genomes:
                continue

            contaminant_genome = self.reference_genomes[species]
            matches = 0

            # Simple k-mer based species identification
            k = 25
            target_kmers = set()
            contaminant_kmers = set()

            for i in range(len(target_genome) - k + 1):
                target_kmers.add(target_genome[i : i + k])

            for i in range(len(contaminant_genome) - k + 1):
                contaminant_kmers.add(contaminant_genome[i : i + k])

            for seq in sequences:
                seq_kmers = set()
                for i in range(len(seq) - k + 1):
                    seq_kmers.add(seq[i : i + k])

                # Check similarity to contaminant vs target
                contaminant_sim = len(seq_kmers.intersection(contaminant_kmers)) / len(seq_kmers) if seq_kmers else 0
                target_sim = len(seq_kmers.intersection(target_kmers)) / len(seq_kmers) if seq_kmers else 0

                if contaminant_sim > target_sim * 1.5:  # Significantly more similar to contaminant
                    matches += 1

            contamination_rate = matches / len(sequences) if sequences else 0

            if contamination_rate > 0.05:  # 5% threshold
                results.append(
                    {
                        "species": species,
                        "contamination_rate": contamination_rate,
                        "matches": matches,
                    }
                )

        detected = len(results) > 0

        return {
            "detected": detected,
            "target_species": target_species,
            "contaminants": results,
            "total_sequences_analyzed": len(sequences),
        }

    def detect_adapter_contamination(
        self, sequences: List[str], adapters: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """Detect adapter sequence contamination.

        Args:
            sequences: List of sequences to analyze
            adapters: List of adapter sequences to check

        Returns:
            Dictionary with adapter contamination results
        """
        if not sequences:
            return {"detected": False, "adapters": []}

        # Default Illumina adapters
        if adapters is None:
            adapters = [
                "AGATCGGAAGAG",  # TruSeq Universal Adapter
                "GATCGGAAGAG",  # TruSeq Adapter, Read 1
                "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",  # TruSeq Adapter, Read 2
                "CTGTCTCTTAT",  # Nextera Transposase Sequence
            ]

        logger.info(f"Detecting adapter contamination in {len(sequences)} sequences")

        results = []
        total_sequences = len(sequences)

        for adapter in adapters:
            adapter_matches = 0
            match_positions = []

            for seq in sequences:
                # Check for exact adapter matches
                if adapter in seq:
                    adapter_matches += 1
                    # Find position
                    pos = seq.find(adapter)
                    match_positions.append(pos)

                # Also check for partial matches at sequence ends
                seq_start = seq[: len(adapter)]
                seq_end = seq[-len(adapter) :]

                if adapter.startswith(seq_start[:10]) or adapter.endswith(seq_end[-10:]):
                    adapter_matches += 1
                    match_positions.append(-1)  # End match

            contamination_rate = adapter_matches / total_sequences

            if contamination_rate > 0.01:  # 1% threshold
                results.append(
                    {
                        "adapter_sequence": adapter,
                        "contamination_rate": contamination_rate,
                        "matches": adapter_matches,
                        "positions": match_positions[:10],  # First 10 positions
                    }
                )

        detected = len(results) > 0

        return {
            "detected": detected,
            "adapters": results,
            "total_sequences_analyzed": total_sequences,
        }

    def detect_duplication_contamination(self, sequences: List[str], max_duplicates: int = 10) -> Dict[str, Any]:
        """Detect excessive sequence duplication that may indicate contamination.

        Args:
            sequences: List of sequences to analyze
            max_duplicates: Maximum allowed duplicates for a sequence

        Returns:
            Dictionary with duplication contamination results
        """
        if not sequences:
            return {"detected": False, "duplicated_sequences": []}

        logger.info(f"Detecting duplication contamination in {len(sequences)} sequences")

        seq_counts = Counter(sequences)
        duplicated_sequences = []

        for seq, count in seq_counts.most_common():
            if count > max_duplicates:
                duplicated_sequences.append(
                    {
                        "sequence": seq[:50] + "..." if len(seq) > 50 else seq,
                        "count": count,
                        "percentage": (count / len(sequences)) * 100,
                    }
                )

        # Calculate duplication metrics
        unique_sequences = len(seq_counts)
        total_sequences = len(sequences)
        duplication_rate = 1 - (unique_sequences / total_sequences)

        detected = len(duplicated_sequences) > 0 or duplication_rate > 0.8

        return {
            "detected": detected,
            "duplication_rate": duplication_rate,
            "unique_sequences": unique_sequences,
            "total_sequences": total_sequences,
            "duplicated_sequences": duplicated_sequences,
            "max_duplicates_threshold": max_duplicates,
        }

    def comprehensive_contamination_analysis(
        self, sequences: List[str], target_species: Optional[str] = None
    ) -> Dict[str, Any]:
        """Perform comprehensive contamination analysis.

        Args:
            sequences: List of sequences to analyze
            target_species: Expected species (for cross-species detection)

        Returns:
            Dictionary with comprehensive contamination analysis results
        """
        logger.info("Performing comprehensive contamination analysis")

        results = {
            "microbial_contamination": self.detect_microbial_contamination(sequences),
            "adapter_contamination": self.detect_adapter_contamination(sequences),
            "duplication_contamination": self.detect_duplication_contamination(sequences),
        }

        if target_species and self.reference_genomes:
            other_species = [s for s in self.reference_genomes.keys() if s != target_species]
            results["cross_species_contamination"] = self.detect_cross_species_contamination(
                sequences, target_species, other_species
            )

        # Overall assessment
        any_detected = any(result.get("detected", False) for result in results.values())

        # Calculate contamination severity score
        severity_scores = []
        for analysis_name, analysis_result in results.items():
            if analysis_result.get("detected", False):
                if "contamination_rate" in analysis_result:
                    severity_scores.append(analysis_result["contamination_rate"] * 100)
                elif "duplication_rate" in analysis_result:
                    severity_scores.append(analysis_result["duplication_rate"] * 100)
                else:
                    severity_scores.append(50)  # Default severity for detected contamination

        overall_severity = np.mean(severity_scores) if severity_scores else 0

        results["summary"] = {
            "contamination_detected": any_detected,
            "overall_severity_score": overall_severity,
            "severity_level": self._classify_severity(overall_severity),
            "analyses_performed": list(results.keys()),
        }

        return results

    def _classify_severity(self, severity_score: float) -> str:
        """Classify contamination severity."""
        if severity_score >= 50:
            return "high"
        elif severity_score >= 20:
            return "moderate"
        elif severity_score >= 5:
            return "low"
        else:
            return "none"




def detect_rna_contamination(dna_sequences: List[str]) -> Dict[str, Any]:
    """Detect RNA contamination in DNA sequencing data."""
    rna_bases = ['U', 'u']
    total_sequences = len(dna_sequences)
    contaminated_sequences = sum(1 for seq in dna_sequences if any(base in seq.upper() for base in rna_bases))
    contamination_rate = contaminated_sequences / total_sequences if total_sequences > 0 else 0
    return {
        'detected': contamination_rate > 0.001,
        'contamination_rate': contamination_rate,
        'contaminated_sequences': contaminated_sequences,
        'total_sequences': total_sequences,
    }


def detect_vector_contamination(sequences: List[str], vector_sequences: Optional[List[str]] = None) -> Dict[str, Any]:
    """Detect vector/plasmid contamination. Returns {seq_idx: [vector_names]}."""
    vector_patterns = {
        'pUC19': 'GCTCTAGAACTAGTGGATC',
        'pBR322': 'GGATCCCCGGGTACCGAGCTCGAATTC',
        'M13': 'GTAAAACGACGGCCAGT',
        'pET': 'TAATACGACTCACTATA',
    }
    results: Dict[str, Any] = {}
    for idx, seq in enumerate(sequences):
        matched = []
        for vector_name, vector_seq in vector_patterns.items():
            if vector_seq in seq.upper() or seq.upper() in vector_seq:
                matched.append(vector_name)
        if not matched and vector_sequences:
            for vs in vector_sequences:
                if vs in seq:
                    matched.append(vs)
        if matched:
            results[str(idx)] = matched
    return results


def detect_adapter_contamination(sequences: List[str], adapters: Optional[List[str]] = None) -> Dict[str, Any]:
    """Detect adapter contamination. Returns {seq_idx: adapter_seq}."""
    if adapters is None:
        adapters = [
            'AGATCGGAAGAG',
            'GATCGGAAGAG',
            'CTGTCTCTTAT',
        ]
    results: Dict[str, Any] = {}
    for idx, seq in enumerate(sequences):
        for adapter in adapters:
            if adapter in seq:
                results[str(idx)] = adapter
                break
    return results


def detect_cross_species_contamination(
    sequences: List[str],
    reference_genomes: Optional[Dict[str, str]] = None,
    threshold: float = 0.5,
    target_species: Optional[str] = None,
    other_species: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Detect cross-species contamination. Returns {seq_idx: species_name}."""
    if not sequences or not reference_genomes:
        return {}
    results: Dict[str, Any] = {}
    for idx, seq in enumerate(sequences):
        best_species = None
        best_score = 0.0
        for species, genome in reference_genomes.items():
            if len(seq) == 0:
                continue
            match_len = 0
            for k in range(len(seq), 0, -1):
                for start in range(len(seq) - k + 1):
                    if seq[start:start + k] in genome:
                        match_len = max(match_len, k)
                        break
                if match_len >= k:
                    break
            score = match_len / len(seq) if len(seq) > 0 else 0
            if score > best_score:
                best_score = score
                best_species = species
        if best_score >= threshold and best_species is not None:
            results[str(idx)] = best_species
    return results


def detect_mycoplasma_contamination(
    sequences: List[str],
    mycoplasma_genome: Optional[str] = None,
) -> Dict[str, Any]:
    """Detect mycoplasma contamination. Returns {seq_idx: True}."""
    patterns = ['TTAAATTTAAATTT', 'AAATTTAAATTTAAATTT']
    results: Dict[str, Any] = {}
    for idx, seq in enumerate(sequences):
        detected = False
        if mycoplasma_genome:
            if seq in mycoplasma_genome or mycoplasma_genome in seq:
                detected = True
            elif len(seq) >= 8:
                for k in range(min(len(seq), len(mycoplasma_genome)), 7, -1):
                    for start in range(len(seq) - k + 1):
                        if seq[start:start + k] in mycoplasma_genome:
                            detected = True
                            break
                    if detected:
                        break
        if not detected:
            for pattern in patterns:
                if pattern in seq:
                    detected = True
                    break
        if detected:
            results[str(idx)] = True
    return results


def detect_rrna_contamination(
    sequences: List[str],
    custom_patterns: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Detect rRNA contamination. Returns {seq_idx: True}."""
    patterns = custom_patterns or [
        'GGAAGGAG',
        'GTGCCAGCAGCCGCGGTAA',
        'GACGGGCGGTGTGT',
        'CCTACGGGAGGCAGCAG',
    ]
    results: Dict[str, Any] = {}
    for idx, seq in enumerate(sequences):
        for pattern in patterns:
            if pattern in seq.upper():
                results[str(idx)] = True
                break
    return results


def generate_contamination_report(
    contamination_results: Dict[str, Any], output_path: Optional[str | Path] = None
) -> str:
    """Generate a contamination analysis report."""
    report_lines = []
    report_lines.append('=' * 60)
    report_lines.append('METAINFORMANT Contamination Analysis Report')
    report_lines.append('=' * 60)
    report_lines.append('')
    all_sample_ids: set = set()
    for category, matches in contamination_results.items():
        if category == 'summary':
            continue
        if isinstance(matches, dict):
            all_sample_ids.update(matches.keys())
    report_lines.append('Summary:')
    report_lines.append(f'Total samples analyzed: {len(all_sample_ids)}')
    report_lines.append('')
    for category, matches in contamination_results.items():
        if category == 'summary':
            continue
        category_label = category.upper()
        report_lines.append(f'{category_label}:')
        if isinstance(matches, dict) and matches:
            report_lines.append(f'  Detected in {len(matches)} sample(s)')
            for sample_id, detail in list(matches.items())[:5]:
                report_lines.append(f'    Sample {sample_id}: {detail}')
        else:
            report_lines.append('  No contamination detected')
        report_lines.append('')
    report = "\n".join(report_lines)
    if output_path:
        output_path = Path(output_path)
        with open(output_path, 'w') as f:
            f.write(report)
        logger.info(f'Contamination report saved to {output_path}')
    return report
