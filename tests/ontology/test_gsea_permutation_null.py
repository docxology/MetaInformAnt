"""Regression tests for the GSEA permutation null in ``metainformant.ontology.pathway_enrichment.enrichment.gsea``.

``gsea`` built its null distribution by shuffling the rank-metric *weights* over the positions of the ranked list
while the gene set stayed where it was. The members of a set concentrated at the top of the ranking still sit at the
top of every permuted list, so the permuted ES is nearly as large as the observed one and the p-value of a strongly
enriched set is close to 1 (about 0.5 for a top-concentrated set and exactly 1.0 for a bottom-concentrated one in the
example below), with an NES stuck at +-1.

``gsea`` now implements the preranked gene-set permutation of Subramanian et al. (2005): the null is the ES of random
gene sets of the same size drawn from the ranked list, the nominal p-value is one-sided in the direction of the
observed ES against the null ES of the same sign (with the +1 correction), and NES is the ES divided by the mean
magnitude of the same-signed null ES.

The ES oracle below is a separate, exact-rational implementation of the Subramanian running sum. The p-value and NES
oracles replay the seeded draws of the random sets and apply the definitions above to the oracle ES values. The
calibration check asks for roughly uniform p-values over random gene sets. No test doubles.
"""

from __future__ import annotations

import math
import random
from fractions import Fraction

import pytest
from scipy import stats

from metainformant.ontology.pathway_enrichment.enrichment import compute_enrichment_score, gsea

N_GENES = 120
NAMES = [f"G{i:03d}" for i in range(N_GENES)]
# a signed rank metric, strongest positive first and strongest negative last
WEIGHTS = [3.0 - 6.0 * i / (N_GENES - 1) for i in range(N_GENES)]
RANKED = list(zip(NAMES, WEIGHTS))


def reference_es(members: set[str]) -> Fraction:
    """Weighted running-sum ES (Subramanian et al. 2005, p = 1) in exact rational arithmetic."""
    in_set = [name in members for name in NAMES]
    n_hit = sum(in_set)
    norm = sum(abs(Fraction(w)) for w, hit in zip(WEIGHTS, in_set) if hit)
    miss_step = Fraction(1, N_GENES - n_hit)
    running, highest, lowest = Fraction(0), Fraction(0), Fraction(0)
    for w, hit in zip(WEIGHTS, in_set):
        running += abs(Fraction(w)) / norm if hit else -miss_step
        highest, lowest = max(highest, running), min(lowest, running)
    return highest if abs(highest) >= abs(lowest) else lowest


def reference_p_and_nes(observed: Fraction, null: list[Fraction]) -> tuple[float, float]:
    """Nominal p-value and NES from the definitions in the module docstring."""
    if observed >= 0:
        same_sign = [e for e in null if e >= 0]
        as_extreme = [e for e in same_sign if e >= observed]
        magnitudes = [e for e in null if e > 0]
    else:
        same_sign = [e for e in null if e < 0]
        as_extreme = [e for e in same_sign if e <= observed]
        magnitudes = [-e for e in null if e < 0]
    p_value = (1 + len(as_extreme)) / (1 + len(same_sign))
    nes = float(observed / (sum(magnitudes) / len(magnitudes))) if magnitudes else math.nan
    return p_value, nes


TOP = sorted(random.Random(1).sample(NAMES[:25], 15))  # concentrated at the top of the ranking
BOTTOM = sorted(random.Random(2).sample(NAMES[-25:], 15))  # concentrated at the bottom
SCATTERED = sorted(random.Random(3).sample(NAMES, 15))
GENE_SETS = {"TOP": TOP, "BOTTOM": BOTTOM, "SCATTERED": SCATTERED}


# ---------------------------------------------------------------------------
# the enrichment score
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("term", sorted(GENE_SETS))
def test_observed_es_matches_the_exact_rational_oracle(term: str) -> None:
    result = {r["term_id"]: r for r in gsea(RANKED, GENE_SETS, n_permutations=9, min_size=5, random_state=0)}[term]

    assert result["es"] == pytest.approx(float(reference_es(set(GENE_SETS[term]))), rel=1e-12)
    assert compute_enrichment_score(NAMES, set(GENE_SETS[term]), weights=WEIGHTS)["es"] == pytest.approx(
        float(reference_es(set(GENE_SETS[term]))), rel=1e-12
    )


def test_es_has_the_expected_sign_for_each_kind_of_set() -> None:
    assert reference_es(set(TOP)) > Fraction(1, 2)
    assert reference_es(set(BOTTOM)) < -Fraction(1, 2)
    assert abs(reference_es(set(SCATTERED))) < Fraction(1, 2)


# ---------------------------------------------------------------------------
# p-values and NES against the gene-set null
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_a_set_concentrated_at_either_end_of_the_ranking_gets_a_small_p_and_a_large_nes(seed: int) -> None:
    results = {r["term_id"]: r for r in gsea(RANKED, GENE_SETS, n_permutations=199, min_size=5, random_state=seed)}

    assert results["TOP"]["p_value"] <= 0.03 and results["TOP"]["nes"] > 2.0
    assert results["BOTTOM"]["p_value"] <= 0.03 and results["BOTTOM"]["nes"] < -2.0
    assert results["SCATTERED"]["p_value"] > 0.1 and abs(results["SCATTERED"]["nes"]) < 2.0
    assert (
        results["TOP"]["fdr"] < results["SCATTERED"]["fdr"] and results["BOTTOM"]["fdr"] < results["SCATTERED"]["fdr"]
    )


def test_p_value_and_nes_replay_exactly_from_the_seeded_random_sets() -> None:
    seed, n_permutations = 5, 60
    results = {
        r["term_id"]: r for r in gsea(RANKED, GENE_SETS, n_permutations=n_permutations, min_size=5, random_state=seed)
    }

    rng = random.Random(seed)  # one generator, consumed set by set in dictionary order, as documented
    for term, genes in GENE_SETS.items():
        size = len(set(genes) & set(NAMES))
        null = [reference_es({NAMES[i] for i in rng.sample(range(N_GENES), size)}) for _ in range(n_permutations)]
        expected_p, expected_nes = reference_p_and_nes(reference_es(set(genes)), null)

        assert results[term]["p_value"] == pytest.approx(expected_p, rel=1e-12), term
        assert results[term]["nes"] == pytest.approx(expected_nes, rel=1e-9), term


def test_p_values_of_random_gene_sets_are_roughly_uniform() -> None:
    set_rng = random.Random(7)
    random_sets = {f"R{j:02d}": set_rng.sample(NAMES, 15) for j in range(80)}

    results = gsea(RANKED, random_sets, n_permutations=99, min_size=5, random_state=11)
    p_values = [r["p_value"] for r in results]

    assert len(p_values) == 80
    assert 0.35 <= sum(p_values) / len(p_values) <= 0.65
    assert stats.kstest(p_values, "uniform").pvalue > 0.01
    assert 0.05 <= sum(p <= 0.2 for p in p_values) / len(p_values) <= 0.40


def test_a_set_that_covers_the_whole_ranked_list_is_never_significant() -> None:
    result = gsea(RANKED, {"ALL": NAMES}, n_permutations=20, min_size=5, random_state=0)[0]

    assert result["es"] == 0.0 and result["p_value"] == 1.0
    assert math.isnan(result["nes"]), "NES is undefined when no permutation produced an ES of the observed sign"


def test_a_set_that_is_missing_from_the_ranking_is_skipped() -> None:
    assert gsea(RANKED, {"ABSENT": ["X1", "X2", "X3", "X4", "X5", "X6"]}, n_permutations=10, min_size=5) == []
