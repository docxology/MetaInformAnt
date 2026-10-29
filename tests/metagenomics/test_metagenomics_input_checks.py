"""Silent fallbacks and placeholder behaviour found while fixing the metagenomics defect list.

Each class pins one behaviour that used to give a plausible but wrong result without saying so: a blank
Krona chart for count data, reverse-strand reads that never classify, output order that follows the string
hash seed, a mismatched error model input that was ignored, and an empty pathway result for a database
that has no definitions. Real implementations.
"""

from __future__ import annotations

import json
import os
import random
import subprocess
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from metainformant.metagenomics.amplicon.asv_denoising import estimate_error_rates
from metainformant.metagenomics.functional.pathways import (
    PathwayDefinition,
    compare_pathway_profiles,
    reconstruct_pathways,
)
from metainformant.metagenomics.shotgun.profiling import KmerIndex, build_kmer_index, profile_community
from metainformant.metagenomics.visualization.plots import plot_krona_chart


def _revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _package_root() -> str:
    import metainformant

    return os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__)))


class TestKronaChartNormalisation:
    """Abundances were normalised twice, so count data drew no wedge at all."""

    FRACTIONS = {
        "Bacteria;Firmicutes;Bacilli": 0.35,
        "Bacteria;Proteobacteria;Gammaproteobacteria": 0.25,
        "Archaea;Euryarchaeota;Methanobacteria": 0.40,
    }

    @staticmethod
    def _wedges(taxonomy: dict[str, float]) -> list[np.ndarray]:
        fig = plot_krona_chart(taxonomy)
        try:
            return [np.asarray(c.get_paths()[0].vertices) for c in fig.axes[0].collections]
        finally:
            plt.close(fig)

    def test_counts_draw_the_same_chart_as_fractions(self) -> None:
        fractions = self._wedges(self.FRACTIONS)
        counts = self._wedges({lineage: value * 1000 for lineage, value in self.FRACTIONS.items()})
        assert len(fractions) == 8  # 2 domains, 3 phyla, 3 classes
        assert len(counts) == len(fractions)
        for a, b in zip(fractions, counts):
            np.testing.assert_allclose(a, b)

    def test_a_child_takes_its_share_of_the_parents_arc(self) -> None:
        """Half of Bacteria stops at the domain, so Firmicutes spans half the arc of its parent, not all of it."""
        wedges = self._wedges({"Bacteria;Firmicutes": 50.0, "Bacteria": 50.0})
        extents = [float(w[:, 0].max() - w[:, 0].min()) for w in wedges]
        assert len(wedges) == 2
        assert extents[0] == pytest.approx(2 * np.pi)  # Bacteria is the whole circle
        assert extents[1] == pytest.approx(np.pi)  # Firmicutes holds half of Bacteria

    def test_small_taxa_are_filtered_by_their_share_of_the_total(self) -> None:
        taxonomy = {"A;a": 990.0, "B;b": 5.0, "C;c": 5.0}  # B and C are 0.5% each: under min_fraction 1%
        assert len(self._wedges(taxonomy)) == 2  # only A and a

    def test_abundances_without_a_positive_total_are_rejected(self) -> None:
        for taxonomy in ({"A;a": 0.0}, {"A;a": 1.0, "B;b": -1.0}):
            with pytest.raises(ValueError, match="positive sum"):
                plot_krona_chart(taxonomy)


class TestProfilingIsStrandIndependent:
    """The index held forward k-mers only, so every read from the other strand went unclassified."""

    @staticmethod
    def _references() -> tuple[dict[str, str], dict[str, list[tuple[str, str]]]]:
        rng = random.Random(14)
        refs = {"ref1": _random_dna(300, rng), "ref2": _random_dna(300, rng)}
        tax = {
            "ref1": [("domain", "Bacteria"), ("phylum", "Firmicutes")],
            "ref2": [("domain", "Bacteria"), ("phylum", "Proteobacteria")],
        }
        return refs, tax

    def test_reads_from_either_strand_are_classified_alike(self) -> None:
        refs, tax = self._references()
        forward = refs["ref1"][40:110]
        for read in (forward, _revcomp(forward)):
            profile = profile_community([read], reference_sequences=refs, reference_taxonomy=tax, k=21)
            assert profile.classified_reads == 1, read == forward
            assert [t.name for t in profile.taxa] == ["Firmicutes"]

    def test_mixed_strands_count_every_read(self) -> None:
        refs, tax = self._references()
        reads = [refs["ref2"][i : i + 60] if i % 20 else _revcomp(refs["ref2"][i : i + 60]) for i in range(0, 200, 10)]
        profile = profile_community(reads, reference_sequences=refs, reference_taxonomy=tax, k=21)
        assert profile.classified_reads == len(reads)
        assert profile.classification_rate == 1.0

    def test_kmer_count_is_the_number_of_distinct_kmers_hit(self) -> None:
        """``TaxonProfile.kmer_count`` was documented but never filled in (always 0)."""
        refs, tax = self._references()
        ref = refs["ref1"]
        # k-mer start positions 10..29 and 20..39 (the second read from the other strand): 30 distinct
        reads = [ref[10:50], _revcomp(ref[20:60])]
        profile = profile_community(reads, reference_sequences=refs, reference_taxonomy=tax, k=21)
        assert [(t.name, t.read_count, t.kmer_count) for t in profile.taxa] == [("Firmicutes", 2, 30)]
        one_read = profile_community(reads[:1], reference_sequences=refs, reference_taxonomy=tax, k=21)
        assert one_read.taxa[0].kmer_count == 20

    def test_index_holds_canonical_kmers(self) -> None:
        refs, tax = self._references()
        index = build_kmer_index(refs, tax, k=21)
        assert index.canonical is True
        expected = {
            min(seq[i : i + 21], _revcomp(seq[i : i + 21])) for seq in refs.values() for i in range(len(seq) - 20)
        }
        assert set(index.kmer_to_taxon) == expected
        assert index.total_kmers == len(expected)

    def test_a_hand_built_forward_index_is_still_queried_as_given(self) -> None:
        sequence = _random_dna(60, random.Random(1))
        kmers = {sequence[i : i + 21]: "t" for i in range(len(sequence) - 20)}
        index = KmerIndex(kmer_to_taxon=kmers, taxon_lineages={"t": [("domain", "X")]}, k=21, total_kmers=len(kmers))
        assert index.canonical is False
        assert profile_community([sequence], database=index).classified_reads == 1


class TestMergedReadsKeepInputOrder:
    """Merged reads were emitted in set order, which follows the per-process string hash seed."""

    def test_output_order_is_the_input_order_for_every_hash_seed(self, tmp_path: Path) -> None:
        rng = random.Random(17)
        forward, reverse = {}, {}
        for i in range(40):
            insert = _random_dna(120, rng)
            forward[f"read_{i:02d}"] = insert[:80]
            reverse[f"read_{i:02d}"] = _revcomp(insert[40:])  # a 40 bp overlap
        code = textwrap.dedent(f"""
            import json, os, metainformant
            assert os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__))) == {_package_root()!r}
            from metainformant.metagenomics.amplicon.asv_denoising import merge_paired_reads
            print(json.dumps(list(merge_paired_reads({forward!r}, {reverse!r}, min_overlap=30))))
            """)
        orders = []
        for hash_seed in ("0", "1", "2"):
            run = subprocess.run(
                [sys.executable, "-c", code],
                env={**os.environ, "PYTHONPATH": _package_root(), "PYTHONHASHSEED": hash_seed},
                cwd=tmp_path,
                capture_output=True,
                text=True,
                check=True,
                timeout=120,
            )
            orders.append(json.loads(run.stdout.strip().splitlines()[-1]))
        assert orders[0] == list(forward)
        assert orders[0] == orders[1] == orders[2]


class TestErrorModelInputs:
    def test_sequences_must_match_the_quality_scores(self) -> None:
        scores = [[30, 30, 30, 30], [30, 30, 30, 20]]
        with pytest.raises(ValueError, match="one entry per read"):
            estimate_error_rates(scores, sequences=["ACGT"])  # used to be ignored without a word
        with pytest.raises(ValueError, match="one entry per read"):
            estimate_error_rates(scores, sequences=["ACGT", "ACGT", "ACGT"])

    def test_matching_sequences_and_no_sequences_still_work(self) -> None:
        scores = [[30, 30, 30, 30], [30, 30, 30, 20]]
        empirical = estimate_error_rates(scores, sequences=["ACGT", "ACGA"])
        uniform = estimate_error_rates(scores)
        assert empirical.positions_modeled == uniform.positions_modeled == 4
        assert estimate_error_rates(scores, sequences=[]).transition_rates == uniform.transition_rates


class TestPathwayDatabase:
    """A database without definitions returned an empty list, indistinguishable from 'no pathway present'."""

    ANNOTATIONS = {"gene1": ["K00844", "K00845"]}

    def test_database_without_definitions_is_an_error(self) -> None:
        with pytest.raises(ValueError, match="No built-in pathway definitions for database 'MetaCyc'"):
            reconstruct_pathways(self.ANNOTATIONS, database="MetaCyc")
        with pytest.raises(ValueError, match="MetaCyc"):
            compare_pathway_profiles({"s": self.ANNOTATIONS}, database="MetaCyc")

    def test_custom_definitions_work_for_any_database_name(self) -> None:
        custom = {"p1": PathwayDefinition("p1", "Custom", database="MetaCyc", required_kos=["K00844", "K00845"])}
        results = reconstruct_pathways(self.ANNOTATIONS, database="MetaCyc", pathway_definitions=custom)
        assert [(r.pathway_id, r.completeness) for r in results] == [("p1", 1.0)]

    def test_kegg_still_has_built_in_definitions(self) -> None:
        results = reconstruct_pathways(self.ANNOTATIONS, database="KEGG")
        assert results and results[0].pathway_id == "map00010"


class TestEstimatorsWithNothingToReport:
    """A taxon that is nowhere, a tie, or an undefined estimator is reported as such, not as a made-up answer."""

    def test_absent_taxon_has_no_associated_group(self) -> None:
        from metainformant.metagenomics.comparative.differential_abundance import indicator_species

        counts = [[10, 0], [12, 0], [0, 0], [1, 0]]
        rows = {r["taxon"]: r for r in indicator_species(counts, [0, 0, 1, 1], ["present", "absent"], 49, seed=1)}
        assert rows["absent"]["associated_group"] is None
        assert rows["absent"]["indicator_value"] == 0.0
        assert rows["present"]["associated_group"] == 0

    def test_equal_group_means_have_no_direction(self) -> None:
        from metainformant.metagenomics.comparative.differential_abundance import effect_size_analysis

        for groups in ([0, 0, 1, 1], [0, 1, 2, 2]):
            rows = effect_size_analysis([[5, 10, 15]] * 4, groups, ["a", "b", "c"])
            assert all(r["direction"] is None for r in rows), groups

    def test_directions_are_still_reported_when_means_differ(self) -> None:
        from metainformant.metagenomics.comparative.differential_abundance import effect_size_analysis

        counts = [[100, 5], [90, 6], [5, 100], [6, 90]]
        rows = {r["taxon"]: r["direction"] for r in effect_size_analysis(counts, [0, 0, 1, 1], ["a", "b"])}
        assert rows == {"a": "enriched_in_0", "b": "enriched_in_1"}

    def test_ace_diverges_when_every_rare_species_is_a_singleton(self) -> None:
        import math

        from metainformant.metagenomics.diversity.metrics import alpha_diversity

        assert alpha_diversity([1, 1, 1, 1], metric="ace")["value"] == math.inf
        assert alpha_diversity([1, 1, 1, 50], metric="ace")["value"] == math.inf  # the abundant species do not help
        finite = alpha_diversity([1, 2, 1, 50], metric="ace")["value"]  # a doubleton gives positive coverage
        assert math.isfinite(finite) and finite > 4
