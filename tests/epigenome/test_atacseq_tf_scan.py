"""Transcription factor motif scan of ATAC-seq peaks.

``find_tf_binding_sites`` used to be simulated: it reported every motif in every peak (``found_in_peak = True``),
invented an expected frequency of 10%, ignored ``genome_fasta`` and returned no coordinates. It now scans the peak
DNA for IUPAC consensus motifs on both strands.

The sequences below are an ``A`` background with motifs planted at known offsets, so every expected site, count and
coordinate is derived by hand in the comments. REAL IMPLEMENTATION.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import pytest

from metainformant.epigenome.assays.atacseq import ATACPeak, find_tf_binding_sites


def _sequence(length: int, plants: List[Tuple[int, str]]) -> str:
    """An ``A`` background of ``length`` bases with each ``(offset, text)`` written over it."""
    bases = ["A"] * length
    for offset, text in plants:
        bases[offset : offset + len(text)] = list(text)
    sequence = "".join(bases)
    assert len(sequence) == length
    return sequence


def _write_fasta(tmp_path: Path, chromosomes: dict[str, str]) -> str:
    fasta = tmp_path / "genome.fa"
    fasta.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in chromosomes.items()))
    return str(fasta)


class TestPlantedMotifs:
    def test_a_motif_is_found_on_both_strands(self) -> None:
        # TTGACA at 40 on the forward strand; TGTCAA at 70 is its reverse complement, i.e. the motif on the minus strand
        peak = ATACPeak("chr1", 1000, 1100, score=7.0)
        sequence = _sequence(100, [(40, "TTGACA"), (70, "TGTCAA")])

        results = find_tf_binding_sites([peak], {"TF": "TTGACA"}, peak_sequences=[sequence])

        sites = results["TF_positions"]
        assert [(s["position"], s["end"], s["strand"]) for s in sites] == [(1040, 1046, "+"), (1070, 1076, "-")]
        # the minus-strand site is read in the motif's orientation
        assert [s["sequence"] for s in sites] == ["TTGACA", "TTGACA"]
        assert all(
            (s["chromosome"], s["peak_start"], s["peak_end"], s["peak_score"]) == ("chr1", 1000, 1100, 7.0)
            for s in sites
        )
        assert results["motif_counts"] == {"TF": 1}  # one peak holds a site
        assert results["site_counts"] == {"TF": 2}

    def test_degenerate_iupac_codes_match_every_base_they_allow(self) -> None:
        # W is A or T: TTGACW matches TTGACA and TTGACT but not TTGACC or TTGACG
        sequence = _sequence(100, [(10, "TTGACA"), (50, "TTGACT"), (80, "TTGACC")])

        results = find_tf_binding_sites([ATACPeak("chr1", 0, 100)], {"TF": "TTGACW"}, peak_sequences=[sequence])

        assert [(s["position"], s["strand"], s["sequence"]) for s in results["TF_positions"]] == [
            (10, "+", "TTGACA"),
            (50, "+", "TTGACT"),
        ]

    def test_only_peaks_that_hold_the_motif_are_counted(self) -> None:
        # The old code counted every peak as a hit
        peaks = [ATACPeak("chr1", 0, 60), ATACPeak("chr1", 500, 560), ATACPeak("chr2", 0, 60)]
        sequences = [_sequence(60, []), _sequence(60, [(20, "GATTACA")]), _sequence(60, [])]

        results = find_tf_binding_sites(peaks, {"TF": "GATTACA"}, peak_sequences=sequences)

        assert results["peaks_analyzed"] == 3
        assert results["peaks_scanned"] == 3
        assert results["motif_counts"] == {"TF": 1}
        assert [(s["chromosome"], s["position"]) for s in results["TF_positions"]] == [("chr1", 520)]

    def test_a_peak_without_any_motif_reports_no_site(self) -> None:
        results = find_tf_binding_sites([ATACPeak("chr1", 0, 50)], {"TF": "GATTACA"}, peak_sequences=["A" * 50])

        assert results["motif_counts"] == {"TF": 0}
        assert results["site_counts"] == {"TF": 0}
        assert results["TF_positions"] == []

    def test_every_motif_gets_its_own_results(self) -> None:
        sequence = _sequence(100, [(10, "GATTACA"), (60, "TTGACA")])

        results = find_tf_binding_sites(
            [ATACPeak("chr1", 0, 100)], {"ONE": "GATTACA", "TWO": "TTGACA"}, peak_sequences=[sequence]
        )

        assert results["motifs_analyzed"] == 2
        assert [s["position"] for s in results["ONE_positions"]] == [10]
        assert [s["position"] for s in results["TWO_positions"]] == [60]
        assert set(results["motif_counts"]) == set(results["motif_enrichment"]) == {"ONE", "TWO"}

    def test_overlapping_occurrences_are_all_reported(self) -> None:
        # AA starts at 0, 1 and 2 in AAAA (TT, its reverse complement, does not occur)
        results = find_tf_binding_sites([ATACPeak("chr1", 0, 4)], {"TF": "AA"}, peak_sequences=["AAAA"])

        assert [(s["position"], s["strand"]) for s in results["TF_positions"]] == [(0, "+"), (1, "+"), (2, "+")]

    def test_every_site_is_kept(self) -> None:
        # The old code kept only the first ten entries per motif
        results = find_tf_binding_sites([ATACPeak("chr1", 0, 40)], {"TF": "AA"}, peak_sequences=["A" * 40])

        assert len(results["TF_positions"]) == 39 == results["site_counts"]["TF"]

    def test_lower_case_sequence_is_scanned(self) -> None:
        results = find_tf_binding_sites([ATACPeak("chr1", 0, 10)], {"TF": "GATC"}, peak_sequences=["aaGaTcaaaa"])

        assert [(s["position"], s["sequence"]) for s in results["TF_positions"]][:1] == [(2, "GATC")]


class TestGenomeFasta:
    def test_sites_are_reported_in_genomic_coordinates(self, tmp_path: Path) -> None:
        # chr1 is 300 bp with TTGACA at 150 and its reverse complement TGTCAA at 220; the peak covers 100-300
        fasta = _write_fasta(tmp_path, {"chr1": _sequence(300, [(150, "TTGACA"), (220, "TGTCAA")])})

        results = find_tf_binding_sites([ATACPeak("chr1", 100, 300)], {"TF": "TTGACA"}, genome_fasta=fasta)

        assert [(s["position"], s["strand"]) for s in results["TF_positions"]] == [(150, "+"), (220, "-")]

    def test_peaks_on_chromosomes_missing_from_the_fasta_are_not_scanned(self, tmp_path: Path) -> None:
        # The old code ignored the FASTA and reported a hit for every peak, including the one on chr9
        fasta = _write_fasta(tmp_path, {"chr1": _sequence(300, [(150, "TTGACA")])})
        peaks = [ATACPeak("chr1", 100, 300), ATACPeak("chr9", 0, 100)]

        results = find_tf_binding_sites(peaks, {"TF": "TTGACA"}, genome_fasta=fasta)

        assert results["peaks_analyzed"] == 2
        assert results["peaks_scanned"] == 1
        assert results["motif_counts"] == {"TF": 1}
        assert all(s["chromosome"] == "chr1" for s in results["TF_positions"])

    def test_a_fasta_without_the_motif_gives_no_hits(self, tmp_path: Path) -> None:
        fasta = _write_fasta(tmp_path, {"chr1": "A" * 300})

        results = find_tf_binding_sites([ATACPeak("chr1", 0, 300)], {"TF": "TTGACA"}, genome_fasta=fasta)

        assert results["motif_counts"] == {"TF": 0}
        assert results["TF_positions"] == []


class TestExpectedSitesAndEnrichment:
    def test_expectation_follows_the_base_composition(self) -> None:
        # ACGTACGT: each base has frequency 1/4. For the motif AC, one window matches on the forward strand with
        # probability 1/16 and on the reverse strand (GT) with 1/16; there are 7 windows, so 7 * 2/16 = 0.875 are
        # expected. AC occurs at 0 and 4, GT (the minus-strand AC) at 2 and 6: 4 observed.
        results = find_tf_binding_sites([ATACPeak("chr1", 0, 8)], {"AC": "AC"}, peak_sequences=["ACGTACGT"])

        assert results["site_counts"] == {"AC": 4}
        assert results["expected_sites"]["AC"] == pytest.approx(0.875)
        assert results["motif_enrichment"]["AC"] == pytest.approx(4 / 0.875)

    def test_palindromic_motif_counts_once_per_strand(self) -> None:
        # ACGT is its own reverse complement: sites at 0 and 4 on each strand, 4 in all. Probability (1/4)^4 per
        # strand and window, 5 windows: 5 * 2/256 expected.
        results = find_tf_binding_sites([ATACPeak("chr1", 0, 8)], {"ACGT": "ACGT"}, peak_sequences=["ACGTACGT"])

        assert [(s["position"], s["strand"]) for s in results["ACGT_positions"]] == [
            (0, "+"),
            (0, "-"),
            (4, "+"),
            (4, "-"),
        ]
        assert results["expected_sites"]["ACGT"] == pytest.approx(5 * 2 / 256)
        assert results["motif_enrichment"]["ACGT"] == pytest.approx(4 / (5 * 2 / 256))

    def test_degenerate_motif_probability_adds_the_allowed_bases(self) -> None:
        # R is A or G: the forward probability of AR over uniform bases is 1/4 * 1/2 = 1/8, and the reverse pattern YT
        # (Y is C or T) has probability 1/8 too: 7 windows * 2/8 = 1.75 expected
        results = find_tf_binding_sites([ATACPeak("chr1", 0, 8)], {"AR": "AR"}, peak_sequences=["ACGTACGT"])

        assert results["expected_sites"]["AR"] == pytest.approx(1.75)

    def test_enrichment_is_unavailable_rather_than_zero(self) -> None:
        # A motif longer than the sequence has no window, and one that needs a base the sequence never contains has an
        # expectation of zero; neither can be divided by
        results = find_tf_binding_sites(
            [ATACPeak("chr1", 0, 5)], {"LONG": "AAAAAAAA", "ABSENT": "GG"}, peak_sequences=["AAAAA"]
        )

        assert results["expected_sites"] == {"LONG": None, "ABSENT": 0.0}
        assert results["motif_enrichment"] == {"LONG": None, "ABSENT": None}

    def test_sequence_with_no_bases_has_no_expectation(self) -> None:
        results = find_tf_binding_sites([ATACPeak("chr1", 0, 6)], {"TF": "AA"}, peak_sequences=["NNNNNN"])

        assert results["expected_sites"] == {"TF": None}
        assert results["motif_enrichment"] == {"TF": None}


class TestInputValidation:
    def test_no_sequence_is_an_error(self) -> None:
        # The old code returned made-up counts here
        with pytest.raises(ValueError, match="sequence"):
            find_tf_binding_sites([ATACPeak("chr1", 0, 100)], {"TF": "TTGACA"})

    def test_no_sequence_is_an_error_even_without_peaks(self) -> None:
        with pytest.raises(ValueError, match="sequence"):
            find_tf_binding_sites([], {"TF": "TTGACA"})

    def test_genome_fasta_and_peak_sequences_are_exclusive(self, tmp_path: Path) -> None:
        fasta = _write_fasta(tmp_path, {"chr1": "A" * 10})
        with pytest.raises(ValueError, match="not both"):
            find_tf_binding_sites(
                [ATACPeak("chr1", 0, 10)], {"TF": "AA"}, genome_fasta=fasta, peak_sequences=["A" * 10]
            )

    def test_one_sequence_per_peak(self) -> None:
        with pytest.raises(ValueError, match="one per peak"):
            find_tf_binding_sites(
                [ATACPeak("chr1", 0, 10), ATACPeak("chr1", 20, 30)], {"TF": "AA"}, peak_sequences=["A" * 10]
            )

    def test_sequence_must_span_the_peak(self) -> None:
        # Otherwise the genomic coordinates of the sites would be wrong
        with pytest.raises(ValueError, match="10 bp long"):
            find_tf_binding_sites([ATACPeak("chr1", 0, 10)], {"TF": "AA"}, peak_sequences=["A" * 9])

    def test_motifs_must_be_iupac(self) -> None:
        with pytest.raises(ValueError, match="IUPAC"):
            find_tf_binding_sites([ATACPeak("chr1", 0, 10)], {"TF": "AAXA"}, peak_sequences=["A" * 10])

    def test_empty_motif_is_rejected(self) -> None:
        with pytest.raises(ValueError, match="empty"):
            find_tf_binding_sites([ATACPeak("chr1", 0, 10)], {"TF": ""}, peak_sequences=["A" * 10])
