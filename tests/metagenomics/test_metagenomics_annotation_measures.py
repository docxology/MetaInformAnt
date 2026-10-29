"""Annotation reports defined measures (profile coverage, consensus identity) and real hit bounds.

``FunctionalAnnotation.confidence`` was ``min(1, best_score / 100)``: the PSSM score is a mean log-odds
per position, so dividing it by 100 measured nothing. The hit also carried placeholder coordinates
(the whole query, target 0-0), and a sequence with nothing to align scored 0.0, which counted as a hit
whenever ``min_score`` was zero or negative. Real implementations; the profiles are built by hand.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

from metainformant.metagenomics.functional.annotation import (
    FunctionalAnnotation,
    annotate_genes,
    classify_gene_families,
)

CONSENSUS = "CDEFGHIK"


def _profile(length: int = 8) -> dict[int, dict[str, float]]:
    """Column p prefers CONSENSUS[p] (score 5.0); alanine is tolerated (score 1.0); anything else scores -1."""
    return {pos: {CONSENSUS[pos]: 5.0, "A": 1.0} for pos in range(length)}


DB = {"KO:K00001": _profile()}


def _annotate(sequence: str, min_score: float = 0.0) -> FunctionalAnnotation:
    return annotate_genes({"g": sequence}, hmm_db=DB, min_score=min_score)[0]


class TestMeasuresAreDefined:
    def test_confidence_is_gone(self) -> None:
        assert "confidence" not in {f.name for f in dataclasses.fields(FunctionalAnnotation)}
        assert not hasattr(_annotate(CONSENSUS), "confidence")

    def test_full_length_exact_match(self) -> None:
        annotation = _annotate(CONSENSUS)
        assert annotation.profile_coverage == 1.0
        assert annotation.consensus_identity == 1.0
        assert annotation.best_hit is not None
        assert annotation.best_hit.score == pytest.approx(5.0)

    def test_substitutions_lower_identity_not_coverage(self) -> None:
        annotation = _annotate("CDXFGHXK")  # 6 of 8 columns carry the consensus residue; X scores -1
        assert annotation.consensus_identity == pytest.approx(6 / 8)
        assert annotation.profile_coverage == 1.0
        assert annotation.best_hit is not None
        assert annotation.best_hit.score == pytest.approx((6 * 5.0 - 2.0) / 8)

    def test_a_tolerated_residue_scores_but_is_not_consensus(self) -> None:
        annotation = _annotate("ADEFGHIK")  # position 0 holds A (score 1.0) where the consensus residue is C
        assert annotation.consensus_identity == pytest.approx(7 / 8)
        assert annotation.best_hit is not None
        assert annotation.best_hit.score == pytest.approx((7 * 5.0 + 1.0) / 8)

    def test_truncated_protein_covers_part_of_the_profile(self) -> None:
        annotation = _annotate("CDEFGH")  # six of eight columns have a residue
        assert annotation.profile_coverage == pytest.approx(6 / 8)
        assert annotation.consensus_identity == 1.0

    def test_measures_follow_the_best_hit_only(self) -> None:
        db = {"KO:K00001": _profile(), "COG:C0001": {pos: {"W": 5.0} for pos in range(8)}}
        annotation = annotate_genes({"g": CONSENSUS}, hmm_db=db, min_score=0.0)[0]
        assert annotation.kegg_orthologs == ["K00001"]
        assert annotation.best_hit is not None and annotation.best_hit.target_id == "KO:K00001"
        assert annotation.consensus_identity == 1.0

    def test_no_hit_means_no_measures(self) -> None:
        annotation = _annotate("WWWWWWWW", min_score=1.0)
        assert annotation.best_hit is None
        assert annotation.profile_coverage is None
        assert annotation.consensus_identity is None
        assert annotate_genes({"g": CONSENSUS}, hmm_db=None)[0].profile_coverage is None

    def test_sparse_profiles_count_defined_positions(self) -> None:
        sparse = {0: {"A": 5.0}, 3: {"D": 5.0}, 9: {"I": 5.0}}  # three defined positions spread over ten
        annotation = annotate_genes({"g": "AWWDWWWWWI"}, hmm_db={"p": sparse}, min_score=0.0)[0]
        assert annotation.profile_coverage == 1.0
        assert annotation.consensus_identity == 1.0
        shorter = annotate_genes({"g": "AWWD"}, hmm_db={"p": sparse}, min_score=0.0)[0]
        assert shorter.profile_coverage == pytest.approx(2 / 3)


class TestHitBoundsAreReal:
    def test_bounds_of_a_match_inside_the_protein(self) -> None:
        hit = _annotate("WW" + CONSENSUS + "WW").best_hit
        assert hit is not None
        assert (hit.query_start, hit.query_end) == (2, 10)  # the protein is 12 long; the old value was (0, 12)
        assert (hit.target_start, hit.target_end) == (0, 8)

    def test_bounds_of_a_truncated_match(self) -> None:
        hit = _annotate("CDEFGH").best_hit
        assert hit is not None
        assert (hit.query_start, hit.query_end) == (0, 6)
        assert (hit.target_start, hit.target_end) == (0, 6)

    def test_scores_stay_uncalibrated(self) -> None:
        hit = _annotate(CONSENSUS).best_hit
        assert hit is not None
        assert math.isnan(hit.e_value) and math.isnan(hit.bias) and math.isnan(hit.domain_e_value)


class TestNothingToAlignIsNotAHit:
    @pytest.mark.parametrize("min_score", [0.0, -5.0])
    def test_empty_protein_is_never_a_hit(self, min_score: float) -> None:
        annotation = _annotate("", min_score=min_score)
        assert annotation.best_hit is None
        assert annotation.kegg_orthologs == []

    def test_empty_profile_is_never_a_hit(self) -> None:
        annotation = annotate_genes({"g": CONSENSUS}, hmm_db={"KO:K00002": {}}, min_score=-5.0)[0]
        assert annotation.best_hit is None

    def test_family_classification_ignores_unalignable_genes(self) -> None:
        result = classify_gene_families(
            {"empty": "", "good": CONSENSUS}, "KEGG", {"KO:K00001": _profile()}, min_score=0.0
        )
        assert result == {"empty": [], "good": ["KO:K00001"]}
