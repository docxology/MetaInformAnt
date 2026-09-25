"""Real-implementation tests for ChIP-seq motif scanning.

Each test builds a tiny deterministic genome FASTA on disk and scans real
ChIPPeak objects against it. Every expected count and coordinate is derived
in the comments. The scan must honestly report zero matches for windows
without a motif and for peaks whose chromosome is missing from the FASTA.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.epigenome.assays.chipseq import ChIPPeak, find_motifs_in_peaks


def _write_genome(tmp_path: Path, chromosomes: dict[str, str]) -> str:
    fasta = tmp_path / "genome.fa"
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in chromosomes.items()))
    return str(fasta)


class TestMotifScanning:
    def test_forward_and_reverse_strand_matches(self, tmp_path: Path):
        # chr1 = "TTGACAGTTTGACA" (14 bp). The motif TTGACA occurs on the
        # forward strand at positions 0 and 8. Its reverse complement
        # (TGTCAA) occurs nowhere, so there are no minus-strand matches.
        fasta = _write_genome(tmp_path, {"chr1": "TTGACAGTTTGACA"})
        peak = ChIPPeak(chromosome="chr1", start=0, end=15, score=10.0)
        results = find_motifs_in_peaks([peak], fasta, ["TTGACA"], window_size=200)

        assert results["total_peaks_analyzed"] == 1
        assert results["motif_counts"]["TTGACA"] == 1  # one peak with matches
        assert results["match_counts"]["TTGACA"] == 2
        entries = results["motif_positions"]["TTGACA"]
        assert [(e["position"], e["strand"]) for e in entries] == [(0, "+"), (8, "+")]
        assert all(e["sequence"] == "TTGACA" for e in entries)
        assert all(
            e["peak_start"] == 0 and e["peak_end"] == 15 and e["peak_score"] == 10.0
            for e in entries
        )

    def test_reverse_strand_only_match(self, tmp_path: Path):
        # chr2 = "AAATGTCAAGGG". The interval 3-8 on the forward strand is
        # TGTCAA, the reverse complement of TTGACA: reading the minus strand
        # there yields the motif. The forward pattern TTGACA occurs nowhere.
        fasta = _write_genome(tmp_path, {"chr2": "AAATGTCAAGGG"})
        peak = ChIPPeak(chromosome="chr2", start=0, end=12, score=5.0)
        results = find_motifs_in_peaks([peak], fasta, ["TTGACA"], window_size=200)

        assert results["motif_counts"]["TTGACA"] == 1
        assert results["match_counts"]["TTGACA"] == 1
        (entry,) = results["motif_positions"]["TTGACA"]
        assert (entry["position"], entry["strand"], entry["sequence"]) == (
            3,
            "-",
            "TTGACA",
        )

    def test_window_without_match_reports_nothing(self, tmp_path: Path):
        # chr3 has no motif-like sequence at all: zero matches, and no
        # fabricated entries for the pattern.
        fasta = _write_genome(tmp_path, {"chr3": "GCGCGCGCGCGC"})
        peak = ChIPPeak(chromosome="chr3", start=0, end=12, score=3.0)
        results = find_motifs_in_peaks([peak], fasta, ["TTGACA"], window_size=200)

        assert results["total_peaks_analyzed"] == 1
        assert "TTGACA" not in results["motif_counts"]
        assert "TTGACA" not in results["match_counts"]
        assert "TTGACA" not in results["motif_positions"]

    def test_peak_on_missing_chromosome_is_skipped(self, tmp_path: Path):
        # chr9 is absent from the FASTA: the peak contributes no matches and
        # none are fabricated. Empty peaks list also works without touching
        # the FASTA.
        fasta = _write_genome(tmp_path, {"chr1": "TTGACAGTTTGACA"})
        peak = ChIPPeak(chromosome="chr9", start=0, end=15, score=1.0)
        results = find_motifs_in_peaks([peak], fasta, ["TTGACA"])

        assert results["total_peaks_analyzed"] == 1
        assert results["motif_counts"] == {}
        assert results["match_counts"] == {}
        assert results["motif_positions"] == {}
        assert find_motifs_in_peaks([], fasta, ["TTGACA"])["total_peaks_analyzed"] == 0

    def test_window_clipped_to_peak_bounds(self, tmp_path: Path):
        # Peak chr1:0-9 keeps only the first 9 bases ("TTGACAGTT"), so only
        # the motif at position 0 is inside the window.
        fasta = _write_genome(tmp_path, {"chr1": "TTGACAGTTTGACA"})
        peak = ChIPPeak(chromosome="chr1", start=0, end=9, score=2.0)
        results = find_motifs_in_peaks([peak], fasta, ["TTGACA"], window_size=200)

        assert results["match_counts"]["TTGACA"] == 1
        (entry,) = results["motif_positions"]["TTGACA"]
        assert (entry["position"], entry["strand"]) == (0, "+")

    def test_iupac_degenerate_code_expansion(self, tmp_path: Path):
        # Pattern TTGACAK (K = G or T): chr1 ends "TTGACA" at 8 with no
        # following base inside the chromosome, but the copy at 0 is
        # followed by G, so exactly one match is reported.
        fasta = _write_genome(tmp_path, {"chr1": "TTGACAGTTTGACA"})
        peak = ChIPPeak(chromosome="chr1", start=0, end=15, score=1.0)
        results = find_motifs_in_peaks([peak], fasta, ["TTGACAK"], window_size=200)

        assert results["match_counts"]["TTGACAK"] == 1
        (entry,) = results["motif_positions"]["TTGACAK"]
        assert (entry["position"], entry["strand"], entry["sequence"]) == (
            0,
            "+",
            "TTGACAG",
        )

    def test_invalid_iupac_pattern_raises(self, tmp_path: Path):
        fasta = _write_genome(tmp_path, {"chr1": "TTGACAGTTTGACA"})
        peak = ChIPPeak(chromosome="chr1", start=0, end=15, score=1.0)
        with pytest.raises(ValueError, match="IUPAC"):
            find_motifs_in_peaks([peak], fasta, ["TTGACAX"])
