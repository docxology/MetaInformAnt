"""Tests for metagenomics amplicon submodule.

Tests OTU clustering, ASV denoising, chimera detection, taxonomy classification.
Uses real implementations -- real-implementation policy.
"""

from __future__ import annotations

import contextlib
import json
import os
import random
import subprocess
import sys
import textwrap
from collections.abc import Iterator

import numpy as np
import pytest
from scipy.stats import poisson

from metainformant.metagenomics.amplicon.asv_denoising import (
    DenoisingResult,
    ErrorModel,
    denoise_sequences,
    estimate_error_rates,
    merge_paired_reads,
)
from metainformant.metagenomics.amplicon.otu_clustering import (
    ClusteringResult,
    calculate_identity,
    cluster_otus,
    filter_chimeras,
)
from metainformant.metagenomics.amplicon.taxonomy import (
    TAXONOMIC_RANKS,
    TaxonomyAssignment,
    TaxonomyNode,
    build_taxonomy_tree,
    calculate_confidence,
    classify_taxonomy,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def identical_sequences() -> dict[str, str]:
    """Three identical sequences."""
    seq = "ATCGATCGATCGATCGATCG"
    return {"s1": seq, "s2": seq, "s3": seq}


@pytest.fixture
def dissimilar_sequences() -> dict[str, str]:
    """Sequences with low identity to each other."""
    return {
        "s1": "ATCGATCGATCGATCGATCG",
        "s2": "TTTTTTTTTTTTTTTTTTTT",
        "s3": "CCCCCCCCCCCCCCCCCCCC",
    }


@pytest.fixture
def reference_taxonomy() -> tuple[dict[str, str], dict[str, list[tuple[str, str]]]]:
    """Reference database with taxonomy."""
    ref_db = {
        "ref1": "ATCGATCGATCGATCGATCGATCG",
        "ref2": "TTTTTTTTGGGGGGGGCCCCCCCC",
        "ref3": "AAAACCCCGGGGTTTTAAAACCCC",
    }
    ref_tax = {
        "ref1": [("domain", "Bacteria"), ("phylum", "Firmicutes"), ("class", "Bacilli")],
        "ref2": [("domain", "Bacteria"), ("phylum", "Proteobacteria"), ("class", "Gammaproteobacteria")],
        "ref3": [("domain", "Archaea"), ("phylum", "Euryarchaeota"), ("class", "Methanobacteria")],
    }
    return ref_db, ref_tax


# ---------------------------------------------------------------------------
# Tests: calculate_identity
# ---------------------------------------------------------------------------


class TestCalculateIdentity:
    """Tests for pairwise sequence identity calculation."""

    def test_identical_sequences_return_one(self) -> None:
        """Identical sequences must have identity = 1.0."""
        assert calculate_identity("ATCGATCG", "ATCGATCG") == 1.0

    def test_different_sequences_less_than_one(self) -> None:
        """Different sequences must have identity < 1.0."""
        identity = calculate_identity("ATCGATCG", "TTTTTTTT")
        assert identity < 1.0

    def test_single_mismatch(self) -> None:
        """One mismatch in 8 bases should give identity near 7/8."""
        identity = calculate_identity("ATCGATCG", "ATCGTTCG")
        assert 0.5 < identity < 1.0

    def test_empty_sequence_raises(self) -> None:
        """Empty sequences should raise ValueError."""
        with pytest.raises(ValueError, match="non-empty"):
            calculate_identity("", "ATCG")

    def test_empty_second_sequence_raises(self) -> None:
        with pytest.raises(ValueError, match="non-empty"):
            calculate_identity("ATCG", "")

    def test_strips_gaps(self) -> None:
        """Gaps should be stripped before comparison."""
        identity = calculate_identity("ATC-GATCG", "ATCGATCG")
        assert identity == 1.0

    def test_case_insensitive(self) -> None:
        """Identity calculation should be case-insensitive."""
        assert calculate_identity("atcgatcg", "ATCGATCG") == 1.0


# ---------------------------------------------------------------------------
# Tests: cluster_otus
# ---------------------------------------------------------------------------


class TestClusterOtus:
    """Tests for OTU clustering."""

    def test_identical_sequences_one_otu(self, identical_sequences: dict[str, str]) -> None:
        """Identical sequences should cluster into one OTU."""
        result = cluster_otus(identical_sequences, threshold=0.97)
        assert isinstance(result, ClusteringResult)
        assert result.num_otus == 1
        assert result.total_sequences == 3

    def test_dissimilar_sequences_multiple_otus(self, dissimilar_sequences: dict[str, str]) -> None:
        """Very dissimilar sequences should form separate OTUs."""
        result = cluster_otus(dissimilar_sequences, threshold=0.97)
        assert result.num_otus > 1
        assert result.total_sequences == 3

    def test_otu_table_populated(self, identical_sequences: dict[str, str]) -> None:
        """OTU table should map centroids to sizes."""
        result = cluster_otus(identical_sequences, threshold=0.97)
        assert len(result.otu_table) == result.num_otus
        assert sum(result.otu_table.values()) == result.total_sequences

    def test_invalid_threshold_raises(self) -> None:
        """Threshold outside (0, 1] should raise ValueError."""
        seqs = {"s1": "ATCGATCG"}
        with pytest.raises(ValueError, match="Threshold"):
            cluster_otus(seqs, threshold=0.0)
        with pytest.raises(ValueError, match="Threshold"):
            cluster_otus(seqs, threshold=1.5)

    def test_empty_sequences_raises(self) -> None:
        """Empty sequences dict should raise ValueError."""
        with pytest.raises(ValueError, match="not be empty"):
            cluster_otus({})

    def test_abundance_sorting(self) -> None:
        """Sorting by abundance should prioritize high-count sequences."""
        seqs = {"rare": "ATCGATCGATCGATCGATCG", "common": "ATCGATCGATCGATCGATCG"}
        abundance = {"rare": 1, "common": 100}
        result = cluster_otus(seqs, threshold=0.97, abundance=abundance, sort_by="abundance")
        assert result.num_otus == 1
        # The centroid should be the common one (sorted first by abundance)
        assert result.otus[0].centroid_id == "common"

    def test_low_threshold_clusters_more(self, dissimilar_sequences: dict[str, str]) -> None:
        """A very low threshold should merge more sequences."""
        high = cluster_otus(dissimilar_sequences, threshold=0.99)
        low = cluster_otus(dissimilar_sequences, threshold=0.50, prefilter=False)
        assert low.num_otus <= high.num_otus


# ---------------------------------------------------------------------------
# Tests: filter_chimeras
# ---------------------------------------------------------------------------


class TestFilterChimeras:
    """Tests for chimera detection."""

    def test_returns_dict_of_bools(self) -> None:
        seqs = {"s1": "AAAAAAAAAA", "s2": "TTTTTTTTTT", "s3": "AAAAATTTTT"}
        result = filter_chimeras(seqs)
        assert isinstance(result, dict)
        assert all(isinstance(v, bool) for v in result.values())
        assert set(result.keys()) == set(seqs.keys())

    def test_most_abundant_not_chimeric(self) -> None:
        """The two most abundant sequences should never be flagged chimeric in de novo mode."""
        seqs = {"s1": "A" * 50, "s2": "T" * 50, "s3": "A" * 25 + "T" * 25}
        abundance = {"s1": 1000, "s2": 500, "s3": 1}
        result = filter_chimeras(seqs, abundance=abundance)
        assert result["s1"] is False
        assert result["s2"] is False

    def test_reference_db_mode(self) -> None:
        """With a reference DB, all sequences are checked against it."""
        seqs = {"q1": "ATCGATCGATCGATCG", "q2": "TTTTTTTTTTTTTTTT"}
        ref_db = {"r1": "ATCGATCGATCGATCG", "r2": "TTTTTTTTTTTTTTTT"}
        result = filter_chimeras(seqs, reference_db=ref_db)
        assert isinstance(result, dict)
        assert len(result) == 2

    def test_empty_sequences_returns_empty(self) -> None:
        assert filter_chimeras({}) == {}


# ---------------------------------------------------------------------------
# Tests: denoise_sequences
# ---------------------------------------------------------------------------


class TestDenoiseSequences:
    """Tests for ASV denoising."""

    def test_basic_denoising(self) -> None:
        reads = {
            "r1": "ATCGATCGATCG",
            "r2": "ATCGATCGATCG",
            "r3": "ATCGATCGATCG",
            "r4": "ATCGATCGATCG",
            "r5": "ATCGATTGATCG",  # one substitution
        }
        result = denoise_sequences(reads)
        assert isinstance(result, DenoisingResult)
        assert result.num_asvs >= 1
        assert result.total_reads >= 4

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            denoise_sequences({})


# ---------------------------------------------------------------------------
# Tests: estimate_error_rates
# ---------------------------------------------------------------------------


class TestEstimateErrorRates:
    """Tests for error rate estimation."""

    def test_basic_operation(self) -> None:
        scores = [[30, 30, 25, 20], [28, 30, 22, 18]]
        model = estimate_error_rates(scores)
        assert isinstance(model, ErrorModel)
        assert len(model.transition_rates) == 4
        assert model.quality_error_rates[30] == pytest.approx(0.001, rel=1e-3)

    def test_with_sequences(self) -> None:
        scores = [[30, 30, 30, 30], [30, 30, 30, 30]]
        seqs = ["ATCG", "ATCG"]
        model = estimate_error_rates(scores, sequences=seqs)
        assert model.positions_modeled == 4

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            estimate_error_rates([])


# ---------------------------------------------------------------------------
# Tests: merge_paired_reads
# ---------------------------------------------------------------------------


class TestMergePairedReads:
    """Tests for paired-end read merging."""

    def test_basic_merge(self) -> None:
        overlap = "ATCGATCGATCG"
        fwd = {"r1": "AAAAAAAAAA" + overlap}
        # reverse complement of overlap region at start of reverse read
        rev = {"r1": overlap + "GGGGGGGGGG"}
        merged = merge_paired_reads(fwd, rev, min_overlap=5)
        # Should be able to merge if the overlap region matches
        assert isinstance(merged, dict)

    def test_no_common_ids_raises(self) -> None:
        with pytest.raises(ValueError, match="No common read IDs"):
            merge_paired_reads({"r1": "ATCG"}, {"r2": "ATCG"})


# ---------------------------------------------------------------------------
# Tests: classify_taxonomy
# ---------------------------------------------------------------------------


class TestClassifyTaxonomy:
    """Tests for taxonomic classification."""

    def test_basic_classification(
        self,
        reference_taxonomy: tuple[dict[str, str], dict[str, list[tuple[str, str]]]],
    ) -> None:
        ref_db, ref_tax = reference_taxonomy
        query = {"q1": "ATCGATCGATCGATCGATCGATCG"}
        results = classify_taxonomy(query, ref_db, ref_tax, method="naive_bayes", bootstrap_n=10)
        assert len(results) == 1
        assert isinstance(results[0], TaxonomyAssignment)
        assert results[0].sequence_id == "q1"

    def test_blast_method(
        self,
        reference_taxonomy: tuple[dict[str, str], dict[str, list[tuple[str, str]]]],
    ) -> None:
        ref_db, ref_tax = reference_taxonomy
        query = {"q1": "ATCGATCGATCGATCGATCGATCG"}
        results = classify_taxonomy(query, ref_db, ref_tax, method="blast")
        assert len(results) == 1
        assert results[0].method == "blast"

    def test_empty_sequences_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            classify_taxonomy({}, {"r1": "ATCG"})

    def test_empty_reference_raises(self) -> None:
        with pytest.raises(ValueError, match="not be empty"):
            classify_taxonomy({"q1": "ATCG"}, {})

    def test_invalid_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown classification"):
            classify_taxonomy({"q1": "ATCG"}, {"r1": "ATCG"}, method="invalid")


# ---------------------------------------------------------------------------
# Tests: build_taxonomy_tree
# ---------------------------------------------------------------------------


class TestBuildTaxonomyTree:
    """Tests for taxonomy tree construction."""

    def test_basic_tree(self) -> None:
        assignments = [
            TaxonomyAssignment("s1", [("domain", "Bacteria"), ("phylum", "Firmicutes")]),
            TaxonomyAssignment("s2", [("domain", "Bacteria"), ("phylum", "Proteobacteria")]),
        ]
        tree = build_taxonomy_tree(assignments)
        assert isinstance(tree, TaxonomyNode)
        assert tree.rank == "root"
        assert tree.total_descendants() == 2
        assert "Bacteria" in tree.children

    def test_tree_to_dict(self) -> None:
        assignments = [
            TaxonomyAssignment("s1", [("domain", "Bacteria")]),
        ]
        tree = build_taxonomy_tree(assignments)
        d = tree.to_dict()
        assert "name" in d
        assert d["total"] == 1


# ---------------------------------------------------------------------------
# Tests: calculate_confidence
# ---------------------------------------------------------------------------


class TestCalculateConfidence:
    """Tests for confidence calculation."""

    def test_returns_float_stats(self) -> None:
        assignments = [
            TaxonomyAssignment("s1", [], confidence={"domain": 0.95, "phylum": 0.80}),
            TaxonomyAssignment("s2", [], confidence={"domain": 0.90, "phylum": 0.60}),
        ]
        stats = calculate_confidence(assignments, min_confidence=0.8)
        assert "domain" in stats
        assert stats["domain"]["mean"] == pytest.approx(0.925, abs=1e-6)
        assert 0.0 <= stats["domain"]["fraction_confident"] <= 1.0

    def test_empty_assignments_returns_empty(self) -> None:
        assert calculate_confidence([]) == {}

    def test_n_is_int(self) -> None:
        """The 'n' statistic should be an integer count of assignments."""
        assignments = [
            TaxonomyAssignment("s1", [], confidence={"domain": 0.95}),
            TaxonomyAssignment("s2", [], confidence={"domain": 0.90}),
        ]
        stats = calculate_confidence(assignments)
        assert stats["domain"]["n"] == 2
        assert isinstance(stats["domain"]["n"], int)


# ---------------------------------------------------------------------------
# Tests: added coverage
# ---------------------------------------------------------------------------


class TestMergePairedReadsSequence:
    """Sequence-level assertions for paired-end merging."""

    def test_merge_produces_expected_sequence(self) -> None:
        """Merged read should be forward flank + overlap + reverse tail."""
        overlap = "ATCGATCGATCG"
        forward_part = "AAAAAAAAAA"
        reverse_part = "GGGGGGGGGG"
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        # Raw reverse read whose reverse complement starts with the overlap
        rev_raw = "".join(complement[b] for b in reversed(overlap + reverse_part))
        fwd = {"r1": forward_part + overlap}
        rev = {"r1": rev_raw}
        merged = merge_paired_reads(fwd, rev, min_overlap=5)
        assert merged == {"r1": forward_part + overlap + reverse_part}


class TestDenoiseIdenticalReads:
    """Denoising of identical reads must collapse to one ASV."""

    def test_identical_reads_single_asv(self) -> None:
        seq = "ATCGATCGATCGATCG"
        reads = {f"r{i}": seq for i in range(4)}
        result = denoise_sequences(reads)
        assert result.num_asvs == 1
        assert result.asvs[0].abundance == 4
        assert result.asvs[0].sequence == seq


class TestClassifyTaxonomyThresholds:
    """Confidence-threshold masking and tree branching behavior."""

    def test_unrelated_query_masked_unclassified(
        self,
        reference_taxonomy: tuple[dict[str, str], dict[str, list[tuple[str, str]]]],
    ) -> None:
        ref_db, ref_tax = reference_taxonomy
        # Poly-T query shares no 11-mers with the ACGT reference sequences
        query = {"q1": "TTTTTTTTTTTTTTTTTTTT"}
        results = classify_taxonomy(query, ref_db, ref_tax, method="blast")
        assert all(name == "unclassified" for _, name in results[0].lineage)

    def test_tree_unclassified_branch_stops_traversal(self) -> None:
        assignments = [
            TaxonomyAssignment("s1", [("domain", "unclassified"), ("phylum", "Firmicutes")]),
        ]
        tree = build_taxonomy_tree(assignments)
        assert "unclassified_domain" in tree.children
        unclassified = tree.children["unclassified_domain"]
        assert unclassified.count == 1
        # Traversal stops at the unclassified node: no deeper ranks recorded
        assert unclassified.children == {}


# ---------------------------------------------------------------------------
# Helpers for the regression tests below
# ---------------------------------------------------------------------------

_COMPLEMENT_BASE = {"A": "C", "C": "G", "G": "T", "T": "A"}


def _random_dna(length: int, rng: random.Random, alphabet: str = "ACGT") -> str:
    return "".join(rng.choice(alphabet) for _ in range(length))


def _substitute(seq: str, positions: "range | list[int]") -> str:
    chars = list(seq)
    for pos in positions:
        chars[pos] = _COMPLEMENT_BASE[chars[pos]]
    return "".join(chars)


def _mutate(seq: str, rng: random.Random, n_sub: int, n_ins: int, n_del: int) -> str:
    chars = list(seq)
    for _ in range(n_sub):
        i = rng.randrange(len(chars))
        chars[i] = rng.choice([c for c in "ACGT" if c != chars[i]])
    for _ in range(n_ins):
        chars.insert(rng.randrange(len(chars) + 1), rng.choice("ACGT"))
    for _ in range(n_del):
        if len(chars) > 1:
            del chars[rng.randrange(len(chars))]
    return "".join(chars)


def _subprocess_env() -> dict[str, str]:
    """Environment for a child interpreter that imports the same ``metainformant`` as this process.

    ``tests/conftest.py`` prepends the checkout's ``src`` to ``PYTHONPATH``, so inheriting the
    variable could silently select a different copy of the package than the one under test.
    """
    import metainformant

    package_root = os.path.dirname(os.path.dirname(os.path.abspath(metainformant.__file__)))
    return {**os.environ, "PYTHONPATH": package_root}


@contextlib.contextmanager
def _environment(**values: str | None) -> Iterator[None]:
    """Temporarily set (or unset, with None) real environment variables."""
    saved = {name: os.environ.get(name) for name in values}
    try:
        for name, value in values.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
        yield
    finally:
        for name, old in saved.items():
            if old is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = old


# ---------------------------------------------------------------------------
# Regression tests: OTU clustering prefilter, threshold resolution, chimeras
# ---------------------------------------------------------------------------


class TestOtuPrefilterIsConservative:
    """The k-mer prefilter used to reject pairs that are above the identity threshold."""

    def test_98_percent_pair_forms_one_otu_with_and_without_prefilter(self) -> None:
        base = _random_dna(400, random.Random(7))
        variant = _substitute(base, range(25, 400, 50))  # 8 substitutions in 400 bp
        assert calculate_identity(base, variant) == pytest.approx(0.98)
        seqs = {"a": base, "b": variant}
        assert cluster_otus(seqs, threshold=0.97, prefilter=False).num_otus == 1
        assert cluster_otus(seqs, threshold=0.97, prefilter=True).num_otus == 1

    def test_prefilter_never_changes_the_clustering(self) -> None:
        rng = random.Random(2024)
        for _ in range(8):
            centers = [_random_dna(rng.choice([50, 70, 90]), rng) for _ in range(rng.randint(2, 4))]
            seqs = {}
            for ci, center in enumerate(centers):
                for j in range(rng.randint(2, 5)):
                    seqs[f"c{ci}_{j}"] = _mutate(
                        center, rng, rng.randint(0, max(1, len(center) // 30)), rng.randint(0, 2), rng.randint(0, 2)
                    )
            for threshold in (0.9, 0.97, 0.99):
                on = cluster_otus(seqs, threshold=threshold, prefilter=True)
                off = cluster_otus(seqs, threshold=threshold, prefilter=False)
                assert [(o.centroid_id, o.member_ids) for o in on.otus] == [
                    (o.centroid_id, o.member_ids) for o in off.otus
                ]

    def test_bound_never_rejects_a_pair_that_reaches_the_threshold(self) -> None:
        from metainformant.metagenomics.amplicon.otu_clustering import _kmer_counts, _prefilter_may_reach_threshold

        rng = random.Random(11)
        qualifying = 0
        for _ in range(200):
            length = rng.choice([20, 40, 80])
            base = _random_dna(length, rng, rng.choice(["ACGT", "AC"]))
            other = _mutate(base, rng, rng.randint(0, max(1, length // 15)), rng.randint(0, 3), rng.randint(0, 3))
            identity = calculate_identity(base, other)
            for threshold in (0.5, 0.8, 0.9, 0.95, 0.97, 0.99, 1.0):
                if identity >= threshold:
                    qualifying += 1
                    assert _prefilter_may_reach_threshold(
                        base, _kmer_counts(base), other, _kmer_counts(other), threshold
                    ), (identity, threshold, length)
        assert qualifying > 200  # the check above is not vacuous

    def test_bound_still_rejects_unrelated_sequences(self) -> None:
        from metainformant.metagenomics.amplicon.otu_clustering import _kmer_counts, _prefilter_may_reach_threshold

        rng = random.Random(5)
        a, b = _random_dna(200, rng), _random_dna(200, rng)
        assert not _prefilter_may_reach_threshold(a, _kmer_counts(a), b, _kmer_counts(b), 0.97)


class TestOtuThresholdResolution:
    """META_OTU_THRESHOLD must not override an explicit argument and must be validated."""

    SIMILAR = {"x": "ATCGATCGATCGATCGATCG", "y": "ATCGATCGATCGATCGATGG"}  # 95% identity

    def test_explicit_threshold_beats_environment(self) -> None:
        with _environment(META_OTU_THRESHOLD="0.5"):
            result = cluster_otus(self.SIMILAR, threshold=0.99)
        assert result.threshold == 0.99
        assert result.num_otus == 2

    def test_environment_is_used_when_no_argument_is_given(self) -> None:
        with _environment(META_OTU_THRESHOLD="0.9"):
            result = cluster_otus(self.SIMILAR)
        assert result.threshold == 0.9
        assert result.num_otus == 1

    def test_default_is_97_percent(self) -> None:
        with _environment(META_OTU_THRESHOLD=None):
            assert cluster_otus(self.SIMILAR).threshold == 0.97

    def test_blank_environment_value_means_unset(self) -> None:
        with _environment(META_OTU_THRESHOLD="  "):
            assert cluster_otus(self.SIMILAR).threshold == 0.97

    @pytest.mark.parametrize("value", ["abc", "7", "0", "-0.2", "nan", "inf", "1.0000001"])
    def test_invalid_environment_value_raises(self, value: str) -> None:
        with _environment(META_OTU_THRESHOLD=value):
            with pytest.raises(ValueError, match="META_OTU_THRESHOLD"):
                cluster_otus(self.SIMILAR)

    def test_explicit_argument_does_not_consult_a_bad_environment(self) -> None:
        with _environment(META_OTU_THRESHOLD="abc"):
            assert cluster_otus(self.SIMILAR, threshold=0.9).threshold == 0.9

    def test_nan_explicit_threshold_rejected(self) -> None:
        with pytest.raises(ValueError, match="Threshold"):
            cluster_otus(self.SIMILAR, threshold=float("nan"))


class TestFilterChimerasMinDivergence:
    """min_divergence (percentage points over the closest parent) used to be ignored."""

    @staticmethod
    def _parents_and_chimera() -> tuple[dict[str, str], str]:
        rng = random.Random(21)
        a, b = _random_dna(100, rng), _random_dna(100, rng)
        return {"A": a, "B": b}, a[:50] + b[50:]

    def test_default_flags_the_chimera(self) -> None:
        refs, chimera = self._parents_and_chimera()
        assert filter_chimeras({"q": chimera}, reference_db=refs) == {"q": True}

    def test_divergence_requirement_is_enforced(self) -> None:
        refs, chimera = self._parents_and_chimera()
        # The model beats the closest parent by 33 identity points
        assert filter_chimeras({"q": chimera}, reference_db=refs, min_divergence=30.0) == {"q": True}
        assert filter_chimeras({"q": chimera}, reference_db=refs, min_divergence=40.0) == {"q": False}

    def test_negative_or_nan_min_divergence_rejected(self) -> None:
        with pytest.raises(ValueError, match="min_divergence"):
            filter_chimeras({"q": "ACGT" * 5}, min_divergence=-1.0)
        with pytest.raises(ValueError, match="min_divergence"):
            filter_chimeras({"q": "ACGT" * 5}, min_divergence=float("nan"))


# ---------------------------------------------------------------------------
# Regression tests: ASV denoising p-value and Phred qualities
# ---------------------------------------------------------------------------


class TestPoissonPValue:
    """The abundance p-value is the Poisson upper tail P(X >= k), via the survival function."""

    @pytest.mark.parametrize("k", [1, 2, 5, 10, 30, 50, 100, 300])
    @pytest.mark.parametrize("lam", [1e-3, 0.05, 0.33, 1.0, 7.0, 40.0, 250.0])
    def test_matches_scipy_survival_function(self, k: int, lam: float) -> None:
        from metainformant.metagenomics.amplicon.asv_denoising import _poisson_pvalue

        assert _poisson_pvalue(k, lam) == pytest.approx(float(poisson.sf(k - 1, lam)), rel=1e-9, abs=0.0)

    def test_tail_far_below_double_precision_rounding_is_kept(self) -> None:
        from metainformant.metagenomics.amplicon.asv_denoising import _poisson_pvalue

        # P(X >= 50) for lambda = 0.33 is about 2e-89; 1 - cdf can only resolve about 1e-16
        assert 1e-90 < _poisson_pvalue(50, 0.33) < 1e-88

    def test_edge_cases(self) -> None:
        from metainformant.metagenomics.amplicon.asv_denoising import _poisson_pvalue

        assert _poisson_pvalue(0, 3.0) == 1.0
        assert _poisson_pvalue(5, 0.0) == 0.0
        assert _poisson_pvalue(0, 0.0) == 1.0


class TestDenoiseKeepsRealVariants:
    """A real variant at 50 reads was merged into a 1000-read parent for many read lengths."""

    @staticmethod
    def _reads(length: int, n_parent: int, n_variant: int) -> dict[str, str]:
        parent = _random_dna(length, random.Random(length))
        variant = _substitute(parent, [length // 2])
        reads = {f"p{i}": parent for i in range(n_parent)}
        reads.update({f"v{i}": variant for i in range(n_variant)})
        return reads

    @pytest.mark.parametrize("length", [24, 40, 250, 300])
    def test_one_snp_variant_at_50_reads_is_its_own_asv(self, length: int) -> None:
        result = denoise_sequences(self._reads(length, 1000, 50))
        assert result.num_asvs == 2
        assert result.reads_removed == 0
        assert sorted(a.abundance for a in result.asvs) == [50, 1000]

    def test_singleton_error_variant_is_still_merged(self) -> None:
        result = denoise_sequences(self._reads(100, 1000, 1))
        assert result.num_asvs == 1
        assert result.reads_removed == 1
        assert result.asvs[0].abundance == 1001


class TestDenoiseUsesPhredScores:
    """Phred scores set the error probability of a mismatch and the ASV quality."""

    LENGTH = 200
    POSITION = 100

    def _reads_and_quals(self, mismatch_q: int, variant_reads: int = 50) -> tuple[dict[str, str], dict[str, list[int]]]:
        parent = _random_dna(self.LENGTH, random.Random(21))
        variant = _substitute(parent, [self.POSITION])
        reads = {f"p{i}": parent for i in range(1000)}
        reads.update({f"v{i}": variant for i in range(variant_reads)})
        quals = {}
        for read_id in reads:
            q = [40] * self.LENGTH
            if read_id.startswith("v"):
                q[self.POSITION] = mismatch_q
            quals[read_id] = q
        return reads, quals

    def test_high_quality_mismatch_is_a_real_variant(self) -> None:
        reads, quals = self._reads_and_quals(mismatch_q=40)
        assert denoise_sequences(reads, quality_scores=quals).num_asvs == 2

    def test_low_quality_mismatch_is_explained_by_errors(self) -> None:
        reads, quals = self._reads_and_quals(mismatch_q=10)
        result = denoise_sequences(reads, quality_scores=quals)
        assert result.num_asvs == 1
        assert result.reads_removed == 50

    def test_asv_quality_is_the_mean_phred_of_its_reads(self) -> None:
        reads, quals = self._reads_and_quals(mismatch_q=10)
        result = denoise_sequences(reads, quality_scores=quals)
        # 1000 parent reads at Q40 everywhere; 50 variant reads at Q40 except one Q10 position
        expected = (1000 * self.LENGTH * 40 + 50 * ((self.LENGTH - 1) * 40 + 10)) / (1050 * self.LENGTH)
        assert result.asvs[0].quality == pytest.approx(expected)

    def test_quality_is_zero_when_no_scores_are_supplied(self) -> None:
        reads, _ = self._reads_and_quals(mismatch_q=10)
        assert all(asv.quality == 0.0 for asv in denoise_sequences(reads).asvs)

    def test_reads_without_scores_fall_back_to_the_model(self) -> None:
        reads, quals = self._reads_and_quals(mismatch_q=10)
        partial = {rid: q for rid, q in quals.items() if rid.startswith("p")}  # no scores for the variant reads
        result = denoise_sequences(reads, quality_scores=partial)
        assert result.num_asvs == 2  # default error model: 50 reads of an error variant are implausible


class TestEstimateErrorRatesUsesPhred:
    """Error probabilities come from the Phred scores, averaged as probabilities."""

    SCORES = [[30, 30, 25, 20], [28, 30, 22, 18]]

    def test_transition_rates_follow_the_mean_error_probability(self) -> None:
        model = estimate_error_rates(self.SCORES)
        flat = [q for row in self.SCORES for q in row]
        mean_error = sum(10 ** (-q / 10) for q in flat) / len(flat)
        for i in range(4):
            assert sum(model.transition_rates[i]) == pytest.approx(1.0)
            assert model.transition_rates[i][i] == pytest.approx(1 - mean_error)
            for j in range(4):
                if i != j:
                    assert model.transition_rates[i][j] == pytest.approx(mean_error / 3)
        # Averaging probabilities is not the same as converting the mean score
        assert mean_error > 10 ** (-(sum(flat) / len(flat)) / 10)

    def test_substitution_shares_come_from_sequences_and_rate_from_phred(self) -> None:
        model = estimate_error_rates([[30, 30, 30, 20], [30, 30, 30, 18]], sequences=["ACGT", "ACGA"])
        eps_t = (10 ** (-20 / 10) + 10 ** (-18 / 10)) / 2
        # consensus T was observed once as T and once as A: all of its error mass goes to A
        t_row = model.transition_rates[3]
        assert t_row[0] == pytest.approx(eps_t)
        assert t_row[1] == 0.0 and t_row[2] == 0.0
        assert t_row[3] == pytest.approx(1 - eps_t)
        # consensus A was never miscalled: Phred Q30 error, spread evenly
        assert model.transition_rates[0][0] == pytest.approx(0.999)
        assert model.transition_rates[0][1] == pytest.approx(0.001 / 3)

    def test_accepts_numpy_arrays(self) -> None:
        from_array = estimate_error_rates(np.array(self.SCORES))
        from_lists = estimate_error_rates(self.SCORES)
        np.testing.assert_allclose(from_array.transition_rates, from_lists.transition_rates)
        assert from_array.positions_modeled == 4


class TestAmpliconWithoutScipy:
    """The package must import and give correct results when SciPy is not installed."""

    CODE = textwrap.dedent("""
        import json, random, sys
        sys.modules["scipy"] = None  # makes `import scipy` raise ImportError

        from metainformant.metagenomics.amplicon import asv_denoising
        assert asv_denoising.HAS_SCIPY is False

        length = 40
        rng = random.Random(length)
        parent = "".join(rng.choice("ACGT") for _ in range(length))
        chars = list(parent)
        chars[length // 2] = {"A": "C", "C": "G", "G": "T", "T": "A"}[chars[length // 2]]
        variant = "".join(chars)
        reads = {f"p{i}": parent for i in range(1000)}
        reads.update({f"v{i}": variant for i in range(50)})
        result = asv_denoising.denoise_sequences(reads)
        print(json.dumps({"num_asvs": result.num_asvs, "p": asv_denoising._poisson_pvalue(50, 0.33)}))
        """)

    def test_exact_poisson_tail_without_scipy(self) -> None:
        run = subprocess.run(
            [sys.executable, "-c", self.CODE],
            env=_subprocess_env(),
            capture_output=True,
            text=True,
            check=True,
            timeout=180,
        )
        out = json.loads(run.stdout.strip().splitlines()[-1])
        assert out["num_asvs"] == 2
        assert out["p"] == pytest.approx(float(poisson.sf(49, 0.33)), rel=1e-9)


# ---------------------------------------------------------------------------
# Regression tests: taxonomy tree counts and naive Bayes abstention
# ---------------------------------------------------------------------------


class TestBuildTaxonomyTreeCountsEachSequenceOnce:
    """Lineages ending in "unclassified" used to be counted twice."""

    def test_total_equals_number_of_assignments(self) -> None:
        assignments = [
            TaxonomyAssignment("s1", [("domain", "Bacteria"), ("phylum", "unclassified")]),
            TaxonomyAssignment("s2", [("domain", "unclassified")]),
            TaxonomyAssignment("s3", [("domain", "Bacteria"), ("phylum", "Firmicutes")]),
        ]
        tree = build_taxonomy_tree(assignments)
        assert tree.total_descendants() == 3
        assert tree.count == 0
        bacteria = tree.children["Bacteria"]
        assert bacteria.count == 0  # s1 is counted in its unclassified branch, not at the parent
        assert bacteria.children["unclassified_phylum"].count == 1
        assert bacteria.children["Firmicutes"].count == 1
        assert tree.children["unclassified_domain"].count == 1

    def test_all_unclassified_counts_each_sequence_once(self) -> None:
        assignments = [
            TaxonomyAssignment(f"s{i}", [(rank, "unclassified") for rank in TAXONOMIC_RANKS]) for i in range(4)
        ]
        tree = build_taxonomy_tree(assignments)
        assert tree.total_descendants() == 4
        assert tree.children["unclassified_domain"].count == 4

    def test_fully_classified_and_truncated_lineages_unchanged(self) -> None:
        assignments = [
            TaxonomyAssignment("s1", [("domain", "Bacteria"), ("phylum", "Firmicutes"), ("class", "Bacilli")]),
            TaxonomyAssignment("s2", [("domain", "Bacteria")]),
            TaxonomyAssignment("s3", []),
        ]
        tree = build_taxonomy_tree(assignments)
        assert tree.total_descendants() == 3
        assert tree.children["Bacteria"].count == 1
        assert tree.children["Bacteria"].children["Firmicutes"].children["Bacilli"].count == 1
        assert tree.count == 1  # the empty lineage stops at the root


class TestNaiveBayesAbstains:
    """RDP-style naive Bayes must abstain without evidence instead of always naming a taxon."""

    @staticmethod
    def _database() -> tuple[dict[str, str], dict[str, list[tuple[str, str]]]]:
        rng = random.Random(5)
        refs = {name: _random_dna(150, rng) for name in ("r1", "r2", "r3", "r4")}
        tax = {
            "r1": [("domain", "Bacteria"), ("phylum", "Firmicutes"), ("class", "Bacilli")],
            "r2": [("domain", "Bacteria"), ("phylum", "Firmicutes"), ("class", "Clostridia")],
            "r3": [("domain", "Bacteria"), ("phylum", "Proteobacteria"), ("class", "Gammaproteobacteria")],
            "r4": [("domain", "Archaea"), ("phylum", "Euryarchaeota"), ("class", "Methanobacteria")],
        }
        return refs, tax

    @staticmethod
    def _disjoint_database() -> tuple[dict[str, str], dict[str, list[tuple[str, str]]]]:
        rng = random.Random(8)
        refs = {"a": _random_dna(150, rng, "AC"), "b": _random_dna(150, rng, "AC")}
        tax = {"a": [("domain", "D1"), ("phylum", "P1")], "b": [("domain", "D1"), ("phylum", "P2")]}
        return refs, tax

    def test_query_without_known_kmers_is_unclassified_with_zero_confidence(self) -> None:
        refs, tax = self._disjoint_database()
        query = _random_dna(100, random.Random(1), "GT")  # shares no 8-mer with the references
        result = classify_taxonomy({"q": query}, refs, tax, method="naive_bayes")[0]
        assert all(name == "unclassified" for _, name in result.lineage)
        assert all(conf == 0.0 for conf in result.confidence.values())

    def test_true_fragments_get_the_reference_lineage(self) -> None:
        refs, tax = self._database()
        for ref_id, seq in refs.items():
            for start in (0, 25, 50):
                result = classify_taxonomy({"q": seq[start : start + 100]}, refs, tax)[0]
                assert result.lineage[:3] == tax[ref_id]
                assert all(result.confidence[rank] >= 0.8 for rank in ("domain", "phylum", "class"))

    def test_unrelated_random_queries_are_not_classified(self) -> None:
        refs, tax = self._database()
        rng = random.Random(77)
        called = 0
        for _ in range(60):
            result = classify_taxonomy({"q": _random_dna(120, rng)}, refs, tax)[0]
            called += result.lineage[0][1] != "unclassified"
        assert called == 0

    def test_single_chance_kmer_is_not_enough_evidence(self) -> None:
        refs, tax = self._disjoint_database()
        base = _random_dna(100, random.Random(2), "GT")
        query = base[:50] + refs["a"][30:38] + base[50:]  # exactly one 8-mer shared with the references
        result = classify_taxonomy({"q": query}, refs, tax)[0]
        assert result.lineage[0] == ("domain", "unclassified")
        assert result.confidence["domain"] < 0.8
        # a lenient threshold reports the same evidence
        lenient = classify_taxonomy({"q": query}, refs, tax, confidence_threshold=0.05)[0]
        assert lenient.lineage[0] == ("domain", "D1")

    def test_lineage_is_cut_at_the_first_uncertain_rank(self) -> None:
        refs, tax = self._database()
        chimera = refs["r1"][:60] + refs["r2"][:60]  # same phylum, different classes
        result = classify_taxonomy({"q": chimera}, refs, tax)[0]
        assert result.lineage[:3] == [("domain", "Bacteria"), ("phylum", "Firmicutes"), ("class", "unclassified")]

    def test_bootstrap_is_reproducible_and_seed_dependent(self) -> None:
        refs, tax = self._disjoint_database()
        base = _random_dna(100, random.Random(2), "GT")
        query = base[:50] + refs["a"][30:38] + base[50:]

        def confidences(seed: int) -> dict[str, float]:
            return classify_taxonomy({"q": query}, refs, tax, seed=seed)[0].confidence

        assert confidences(3) == confidences(3)
        assert len({tuple(confidences(seed).values()) for seed in range(6)}) > 1

    def test_classification_does_not_consume_global_random_state(self) -> None:
        refs, tax = self._database()
        random.seed(123)
        expected = random.random()
        random.seed(123)
        classify_taxonomy({"q": refs["r1"][:100]}, refs, tax)
        assert random.random() == expected

    def test_results_do_not_depend_on_other_queries(self) -> None:
        refs, tax = self._disjoint_database()
        base = _random_dna(100, random.Random(2), "GT")
        weak = base[:50] + refs["a"][30:38] + base[50:]
        alone = classify_taxonomy({"w": weak}, refs, tax)[0]
        batch = classify_taxonomy({"x": base, "w": weak}, refs, tax)[1]
        assert alone.confidence == batch.confidence
