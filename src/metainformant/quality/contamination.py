"""Contamination detection and analysis.

This module provides algorithms for detecting various types of contamination
in biological sequencing data, including microbial contamination, cross-species
contamination, and technical artifacts.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Set
import numpy as np
from collections import Counter, defaultdict

from metainformant.core import logging, errors, validation

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
                contaminant_kmers.add(contaminant_seq[i:i+k])

            for seq in sequences:
                seq_kmers = set()
                for i in range(len(seq) - k + 1):
                    seq_kmers.add(seq[i:i+k])

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
                results.append({
                    "contaminant": contaminant_name,
                    "contamination_rate": contamination_rate,
                    "matches": matches,
                    "average_similarity": np.mean(match_positions) if match_positions else 0,
                })

        detected = len(results) > 0

        return {
            "detected": detected,
            "contaminants": results,
            "total_sequences_analyzed": total_sequences,
            "threshold": threshold,
        }

    def detect_cross_species_contamination(self, sequences: List[str],
                                        target_species: str,
                                        other_species: List[str]) -> Dict[str, Any]:
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
                target_kmers.add(target_genome[i:i+k])

            for i in range(len(contaminant_genome) - k + 1):
                contaminant_kmers.add(contaminant_genome[i:i+k])

            for seq in sequences:
                seq_kmers = set()
                for i in range(len(seq) - k + 1):
                    seq_kmers.add(seq[i:i+k])

                # Check similarity to contaminant vs target
                contaminant_sim = len(seq_kmers.intersection(contaminant_kmers)) / len(seq_kmers) if seq_kmers else 0
                target_sim = len(seq_kmers.intersection(target_kmers)) / len(seq_kmers) if seq_kmers else 0

                if contaminant_sim > target_sim * 1.5:  # Significantly more similar to contaminant
                    matches += 1

            contamination_rate = matches / len(sequences) if sequences else 0

            if contamination_rate > 0.05:  # 5% threshold
                results.append({
                    "species": species,
                    "contamination_rate": contamination_rate,
                    "matches": matches,
                })

        detected = len(results) > 0

        return {
            "detected": detected,
            "target_species": target_species,
            "contaminants": results,
            "total_sequences_analyzed": len(sequences),
        }

    def detect_adapter_contamination(self, sequences: List[str],
                                   adapters: Optional[List[str]] = None) -> Dict[str, Any]:
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
                "GATCGGAAGAG",   # TruSeq Adapter, Read 1
                "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",  # TruSeq Adapter, Read 2
                "CTGTCTCTTAT",   # Nextera Transposase Sequence
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
                seq_start = seq[:len(adapter)]
                seq_end = seq[-len(adapter):]

                if adapter.startswith(seq_start[:10]) or adapter.endswith(seq_end[-10:]):
                    adapter_matches += 1
                    match_positions.append(-1)  # End match

            contamination_rate = adapter_matches / total_sequences

            if contamination_rate > 0.01:  # 1% threshold
                results.append({
                    "adapter_sequence": adapter,
                    "contamination_rate": contamination_rate,
                    "matches": adapter_matches,
                    "positions": match_positions[:10],  # First 10 positions
                })

        detected = len(results) > 0

        return {
            "detected": detected,
            "adapters": results,
            "total_sequences_analyzed": total_sequences,
        }

    def detect_duplication_contamination(self, sequences: List[str],
                                       max_duplicates: int = 10) -> Dict[str, Any]:
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
                duplicated_sequences.append({
                    "sequence": seq[:50] + "..." if len(seq) > 50 else seq,
                    "count": count,
                    "percentage": (count / len(sequences)) * 100,
                })

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

    def comprehensive_contamination_analysis(self, sequences: List[str],
                                          target_species: Optional[str] = None) -> Dict[str, Any]:
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
    """Detect RNA contamination in DNA sequencing data.

    Args:
        dna_sequences: List of DNA sequences that should not contain RNA

    Returns:
        Dictionary with RNA contamination detection results
    """
    rna_bases = ['U', 'u']
    total_sequences = len(dna_sequences)
    contaminated_sequences = 0

    for seq in dna_sequences:
        if any(base in seq.upper() for base in rna_bases):
            contaminated_sequences += 1

    contamination_rate = contaminated_sequences / total_sequences if total_sequences > 0 else 0

    return {
        "detected": contamination_rate > 0.001,  # 0.1% threshold
        "contamination_rate": contamination_rate,
        "contaminated_sequences": contaminated_sequences,
        "total_sequences": total_sequences,
    }


def detect_vector_contamination(sequences: List[str],
                              vector_sequences: Optional[List[str]] = None) -> Dict[str, Any]:
    """Detect vector/plasmid contamination.

    Args:
        sequences: List of sequences to analyze
        vector_sequences: List of known vector sequences

    Returns:
        Dictionary with vector contamination detection results
    """
    if not vector_sequences:
        # Common vector elements
        vector_sequences = [
            "GGATCC",  # BamHI site
            "GAATTC",  # EcoRI site
            "CTCGAG",  # XhoI site
            "GTCGAC",  # SalI site
            "CCCGGG",  # SmaI site
        ]

    total_sequences = len(sequences)
    vector_matches = 0
    matched_vectors = defaultdict(int)

    for seq in sequences:
        for vector_seq in vector_sequences:
            if vector_seq in seq:
                vector_matches += 1
                matched_vectors[vector_seq] += 1
                break  # Count each sequence only once

    contamination_rate = vector_matches / total_sequences if total_sequences > 0 else 0

    return {
        "detected": contamination_rate > 0.01,  # 1% threshold
        "contamination_rate": contamination_rate,
        "vector_matches": vector_matches,
        "total_sequences": total_sequences,
        "matched_vectors": dict(matched_vectors),
    }


def generate_contamination_report(contamination_results: Dict[str, Any],
                                output_path: Optional[str | Path] = None) -> str:
    """Generate a contamination analysis report.

    Args:
        contamination_results: Results from contamination analysis
        output_path: Optional path to save the report

    Returns:
        Formatted contamination report as string
    """
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("CONTAMINATION ANALYSIS REPORT")
    report_lines.append("=" * 60)
    report_lines.append("")

    summary = contamination_results.get("summary", {})
    if summary:
        report_lines.append(f"Overall Contamination Detected: {summary.get('contamination_detected', False)}")
        report_lines.append(f"Severity Level: {summary.get('severity_level', 'unknown')}")
        report_lines.append(f"Severity Score: {summary.get('overall_severity_score', 0):.1f}")
        report_lines.append("")

    # Detailed results
    for analysis_name, analysis_result in contamination_results.items():
        if analysis_name == "summary":
            continue

        report_lines.append(f"{analysis_name.replace('_', ' ').title()}:")
        detected = analysis_result.get("detected", False)
        report_lines.append(f"  Detected: {detected}")

        if detected:
            if "contamination_rate" in analysis_result:
                rate = analysis_result["contamination_rate"]
                report_lines.append(f"  Contamination Rate: {rate:.2%}")

            if "contaminants" in analysis_result:
                contaminants = analysis_result["contaminants"]
                report_lines.append(f"  Found {len(contaminants)} contaminant(s):")
                for contaminant in contaminants[:5]:  # Show top 5
                    report_lines.append(f"    - {contaminant}")

            if "adapters" in analysis_result:
                adapters = analysis_result["adapters"]
                report_lines.append(f"  Found {len(adapters)} adapter(s):")
                for adapter in adapters[:3]:  # Show top 3
                    seq = adapter.get("adapter_sequence", "")[:20] + "..."
                    rate = adapter.get("contamination_rate", 0)
                    report_lines.append(f"    - {seq}: {rate:.2%}")

        report_lines.append("")

    # Recommendations
    report_lines.append("Recommendations:")
    severity = summary.get("severity_level", "none")

    if severity == "high":
        report_lines.append("  ✗ High contamination detected. Data quality is severely compromised.")
        report_lines.append("    - Recommend re-sequencing with improved protocols")
        report_lines.append("    - Consider decontamination procedures")
    elif severity == "moderate":
        report_lines.append("  ⚠ Moderate contamination detected.")
        report_lines.append("    - Review sequencing protocols and laboratory procedures")
        report_lines.append("    - Consider filtering out contaminated reads")
    elif severity == "low":
        report_lines.append("  ✓ Low level contamination detected.")
        report_lines.append("    - Monitor in downstream analysis")
        report_lines.append("    - Generally acceptable for most analyses")
    else:
        report_lines.append("  ✓ No significant contamination detected.")
        report_lines.append("    - Data quality appears good")

    report = "\n".join(report_lines)

    if output_path:
        output_path = Path(output_path)
        with open(output_path, 'w') as f:
            f.write(report)
        logger.info(f"Contamination report saved to {output_path}")

    return report
