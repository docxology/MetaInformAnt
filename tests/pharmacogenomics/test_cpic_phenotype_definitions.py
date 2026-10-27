"""CPIC phenotype definitions and the single source of truth for allele values and phenotypes.

Every expected phenotype here is written from the CPIC publication named next to it, not
derived from the package's own tables, so these tests fail if a phenotype band, an allele
value or a second copy of either drifts away from CPIC. A second oracle is built from the
allele *function labels* in ``star_allele`` ("No function", "Decreased function", ...),
which are independent of the numeric activity table the phenotype code reads.

REAL IMPLEMENTATION -- all tests use real implementations.
"""

from __future__ import annotations

import itertools

import pytest

from metainformant.pharmacogenomics.alleles.diplotype import (
    _ACTIVITY_SCORE_TABLES,
    determine_diplotype,
)
from metainformant.pharmacogenomics.alleles.phenotype import (
    _PHENOTYPE_THRESHOLDS,
    MetabolizerPhenotype,
    classify_phenotype,
    get_phenotype_thresholds,
    population_phenotype_frequencies,
    predict_metabolizer_status,
)
from metainformant.pharmacogenomics.alleles.star_allele import _BUILTIN_ALLELE_DEFINITIONS
from metainformant.pharmacogenomics.metabolism.metabolizer_status import (
    classify_metabolizer,
    default_allele_function_table,
    dose_adjustment,
    predict_metabolizer_status as predict_metabolizer,
)

PM = MetabolizerPhenotype.POOR
IM = MetabolizerPhenotype.INTERMEDIATE
NM = MetabolizerPhenotype.NORMAL
RM = MetabolizerPhenotype.RAPID
UM = MetabolizerPhenotype.ULTRARAPID

# (gene, diplotype, expected phenotype), grouped by the CPIC table the expectation comes from.
CPIC_DIPLOTYPES = [
    # CYP2D6 -- CPIC/DPWG 2019 consensus (Caudle et al., Clin Transl Sci 2020): activity score
    # 0 PM; >0 to <1.25 IM; 1.25 to 2.25 NM; >2.25 UM. *10 counts 0.25, *41 0.5.
    ("CYP2D6", "*1/*1", NM),
    ("CYP2D6", "*1/*2", NM),
    ("CYP2D6", "*1/*4", IM),  # AS 1.0
    ("CYP2D6", "*1/*5", IM),  # AS 1.0
    ("CYP2D6", "*41/*41", IM),  # AS 1.0
    ("CYP2D6", "*1/*10", NM),  # AS 1.25
    ("CYP2D6", "*1/*41", NM),  # AS 1.5
    ("CYP2D6", "*4/*41", IM),  # AS 0.5
    ("CYP2D6", "*10/*10", IM),  # AS 0.5
    ("CYP2D6", "*4/*10", IM),  # AS 0.25
    ("CYP2D6", "*4/*4", PM),
    ("CYP2D6", "*5/*5", PM),
    ("CYP2D6", "*1/*1xN", UM),  # one extra copy of *1: AS 3.0
    ("CYP2D6", "*4/*1xN", NM),  # AS 2.0
    # CYP2C19 -- Moriyama et al., Clin Pharmacol Ther 2017 Table 1 and Lee et al., Clin Pharmacol
    # Ther 2022 Table 1: by allele function.
    ("CYP2C19", "*17/*17", UM),
    ("CYP2C19", "*1/*17", RM),
    ("CYP2C19", "*1/*1", NM),
    ("CYP2C19", "*1/*2", IM),
    ("CYP2C19", "*1/*3", IM),
    ("CYP2C19", "*2/*17", IM),
    ("CYP2C19", "*2/*2", PM),
    ("CYP2C19", "*2/*3", PM),
    ("CYP2C19", "*3/*3", PM),
    # CYP2C19 "likely" categories (decreased-function *9/*10), reported here as IM / PM
    ("CYP2C19", "*1/*9", IM),
    ("CYP2C19", "*9/*17", IM),  # sums to 2.0, exactly like *1/*1
    ("CYP2C19", "*9/*9", PM),  # sums to 1.0, exactly like *1/*2
    ("CYP2C19", "*2/*9", PM),
    # SLCO1B1 -- Cooper-DeHoff et al., Clin Pharmacol Ther 2022: normal / decreased / poor function.
    ("SLCO1B1", "*1/*1", NM),
    ("SLCO1B1", "*1/*1B", NM),
    ("SLCO1B1", "*1B/*1B", NM),
    ("SLCO1B1", "*1/*5", IM),
    ("SLCO1B1", "*1/*15", IM),
    ("SLCO1B1", "*1B/*5", IM),
    ("SLCO1B1", "*5/*5", PM),
    ("SLCO1B1", "*5/*15", PM),
    ("SLCO1B1", "*15/*15", PM),
    # DPYD -- Amstutz et al., Clin Pharmacol Ther 2018: activity score 2 NM; 1-1.5 IM; 0-0.5 PM.
    ("DPYD", "*1/*1", NM),
    ("DPYD", "*1/*2A", IM),
    ("DPYD", "*1/*13", IM),
    ("DPYD", "*1/c.2846A>T", IM),  # AS 1.5
    ("DPYD", "*1/HapB3", IM),  # AS 1.5
    ("DPYD", "c.2846A>T/HapB3", IM),  # AS 1.0
    ("DPYD", "*2A/c.2846A>T", PM),  # AS 0.5
    ("DPYD", "*2A/HapB3", PM),  # AS 0.5
    ("DPYD", "*2A/*2A", PM),
    ("DPYD", "*2A/*13", PM),
    # CYP2C9 -- Karnes et al., Clin Pharmacol Ther 2021 and Theken et al. 2020: AS 2 NM; 1-1.5 IM; 0-0.5 PM.
    ("CYP2C9", "*1/*1", NM),
    ("CYP2C9", "*1/*2", IM),  # AS 1.5
    ("CYP2C9", "*1/*11", IM),  # AS 1.5
    ("CYP2C9", "*1/*3", IM),  # AS 1.0
    ("CYP2C9", "*2/*2", IM),  # AS 1.0
    ("CYP2C9", "*2/*3", PM),  # AS 0.5
    ("CYP2C9", "*3/*3", PM),
    # CYP3A5 -- Birdwell et al., Clin Pharmacol Ther 2015: *1/*1 NM; *1 with one no-function allele IM; two PM.
    ("CYP3A5", "*1/*1", NM),
    ("CYP3A5", "*1/*3", IM),
    ("CYP3A5", "*1/*6", IM),
    ("CYP3A5", "*1/*7", IM),
    ("CYP3A5", "*3/*3", PM),
    ("CYP3A5", "*3/*6", PM),
    # TPMT / NUDT15 -- Relling et al., Clin Pharmacol Ther 2019.
    ("TPMT", "*1/*1", NM),
    ("TPMT", "*1/*2", IM),
    ("TPMT", "*1/*3A", IM),
    ("TPMT", "*1/*3C", IM),
    ("TPMT", "*3A/*3A", PM),
    ("TPMT", "*3A/*3C", PM),
    ("TPMT", "*2/*3A", PM),
    ("NUDT15", "*1/*1", NM),
    ("NUDT15", "*1/*2", IM),
    ("NUDT15", "*1/*3", IM),
    ("NUDT15", "*3/*3", PM),
    ("NUDT15", "*2/*3", PM),
    # UGT1A1 -- Gammal et al., Clin Pharmacol Ther 2016 Table 2.
    ("UGT1A1", "*1/*1", NM),
    ("UGT1A1", "*1/*36", NM),
    ("UGT1A1", "*1/*28", IM),
    ("UGT1A1", "*1/*37", IM),
    ("UGT1A1", "*28/*28", PM),
    ("UGT1A1", "*28/*37", PM),
    ("UGT1A1", "*37/*37", PM),
]

_IDS = [f"{gene}-{diplotype}" for gene, diplotype, _ in CPIC_DIPLOTYPES]


def _token(phenotype: MetabolizerPhenotype) -> str:
    return phenotype.value.split()[0].lower()


class TestKnownDiplotypes:
    """Diplotype -> phenotype pairs published by CPIC."""

    @pytest.mark.parametrize(("gene", "diplotype", "expected"), CPIC_DIPLOTYPES, ids=_IDS)
    def test_classify_phenotype_matches_cpic(self, gene: str, diplotype: str, expected: MetabolizerPhenotype) -> None:
        result = classify_phenotype(diplotype, gene)
        assert result["phenotype"] == expected
        assert result["phenotype_abbreviation"] == expected.abbreviation

    @pytest.mark.parametrize(("gene", "diplotype", "expected"), CPIC_DIPLOTYPES, ids=_IDS)
    def test_predict_metabolizer_matches_cpic(self, gene: str, diplotype: str, expected: MetabolizerPhenotype) -> None:
        result = predict_metabolizer({"diplotype": diplotype}, gene)
        assert result["phenotype"] == _token(expected)

    def test_cyp2d6_star1_star4_is_intermediate_in_both_apis(self) -> None:
        # The original defect: classify_phenotype said Normal, predict_metabolizer said intermediate.
        assert classify_phenotype("*1/*4", "CYP2D6")["phenotype"] == IM
        assert predict_metabolizer({"diplotype": "*1/*4"}, "CYP2D6")["phenotype"] == "intermediate"

    def test_cyp2c19_star1_star1_is_normal_not_rapid(self) -> None:
        assert classify_phenotype("*1/*1", "CYP2C19")["phenotype"] == NM

    def test_cyp2c19_star1_star17_is_rapid_not_ultrarapid(self) -> None:
        assert classify_phenotype("*1/*17", "CYP2C19")["phenotype"] == RM

    def test_heterozygous_dpyd_decreased_function_is_not_normal(self) -> None:
        # Treated as Normal, a *1/c.2846A>T carrier would get full-dose fluoropyrimidine.
        assert classify_phenotype("*1/c.2846A>T", "DPYD")["phenotype"] == IM


# CPIC allele functionality records (api.cpicpgx.org/v1/allele, retrieved 2026-10-27) for
# the alleles whose values this package once recorded differently.
CPIC_ACTIVITY_VALUES = {
    ("CYP2D6", "*9"): 0.25,
    ("CYP2D6", "*10"): 0.25,
    ("CYP2D6", "*14"): 0.5,
    ("CYP2D6", "*41"): 0.25,
    ("CYP2C9", "*5"): 0.5,
    ("CYP2C9", "*8"): 0.5,
}
CPIC_FUNCTION_CLASSES = {
    ("SLCO1B1", "*5"): "No function",
    ("SLCO1B1", "*15"): "No function",
    ("UGT1A1", "*37"): "Decreased function",
}


class TestCpicAlleleRecords:
    """Allele values and function classes match the CPIC allele functionality table."""

    @pytest.mark.parametrize(("gene", "allele"), sorted(CPIC_ACTIVITY_VALUES))
    def test_activity_values(self, gene: str, allele: str) -> None:
        assert _ACTIVITY_SCORE_TABLES[gene][allele] == CPIC_ACTIVITY_VALUES[(gene, allele)]
        definition = _BUILTIN_ALLELE_DEFINITIONS[gene].get(allele)
        if definition is not None:
            assert definition["activity_value"] == CPIC_ACTIVITY_VALUES[(gene, allele)]

    @pytest.mark.parametrize(("gene", "allele"), sorted(CPIC_FUNCTION_CLASSES))
    def test_function_classes(self, gene: str, allele: str) -> None:
        value = {"No function": 0.0, "Decreased function": 0.5}[CPIC_FUNCTION_CLASSES[(gene, allele)]]
        assert _ACTIVITY_SCORE_TABLES[gene][allele] == value
        definition = _BUILTIN_ALLELE_DEFINITIONS.get(gene, {}).get(allele)
        if definition is not None:
            assert definition["function"] == CPIC_FUNCTION_CLASSES[(gene, allele)]

    @pytest.mark.parametrize(
        ("gene", "diplotype", "score", "phenotype"),
        [
            ("CYP2D6", "*1/*14", 1.5, NM),
            ("CYP2D6", "*4/*41", 0.25, IM),
            ("CYP2D6", "*1/*41", 1.25, NM),
            ("CYP2C9", "*5/*5", 1.0, IM),
            ("CYP2C9", "*2/*8", 1.0, IM),
            ("CYP2C9", "*3/*5", 0.5, PM),
            ("SLCO1B1", "*1/*5", 1.0, IM),
            ("SLCO1B1", "*5/*15", 0.0, PM),
            ("UGT1A1", "*1/*37", 1.5, IM),
        ],
    )
    def test_scores_and_phenotypes(
        self, gene: str, diplotype: str, score: float, phenotype: MetabolizerPhenotype
    ) -> None:
        result = classify_phenotype(diplotype, gene)
        assert result["activity_score"] == score
        assert result["phenotype"] == phenotype
        assert predict_metabolizer_status(score, gene) == phenotype


class TestAlleleFunctionOracle:
    """Check every diplotype of every gene that has calling definitions against CPIC's rules.

    The rules are applied to the function label (or activity value) of each allele in
    ``star_allele._BUILTIN_ALLELE_DEFINITIONS``, which the phenotype code does not read.
    """

    @staticmethod
    def _classes(gene: str) -> dict[str, str]:
        label_to_class = {
            "Normal function": "normal",
            "Decreased function": "decreased",
            "No function": "none",
            "Increased function": "increased",
        }
        return {name: label_to_class[d["function"]] for name, d in _BUILTIN_ALLELE_DEFINITIONS[gene].items()}

    @staticmethod
    def _by_function(gene: str, first: str, second: str) -> MetabolizerPhenotype:
        pair = [first, second]
        none, decreased = pair.count("none"), pair.count("decreased")
        normal, increased = pair.count("normal"), pair.count("increased")
        if gene in {"TPMT", "NUDT15", "CYP3A5"}:
            return {0: NM, 1: IM, 2: PM}[none]
        if gene == "SLCO1B1":
            # Cooper-DeHoff 2022: one reduced (no or decreased function) allele beside a
            # normal one is decreased function, two are poor function.
            return {0: NM, 1: IM, 2: PM}[none + decreased]
        assert gene == "CYP2C19"
        if increased == 2:
            return UM
        if normal == 1 and increased == 1:
            return RM
        if normal == 2:
            return NM
        if none == 2:
            return PM
        if none == 1 and (normal == 1 or increased == 1):
            return IM
        if decreased == 1 and (normal == 1 or increased == 1):
            return IM  # CPIC "likely intermediate"
        assert decreased == 2 or (decreased == 1 and none == 1)
        return PM  # CPIC "likely poor"

    @staticmethod
    def _by_activity_score(gene: str, total: float) -> MetabolizerPhenotype:
        if gene == "CYP2D6":
            if total == 0:
                return PM
            if total < 1.25:
                return IM
            return NM if total <= 2.25 else UM
        assert gene in {"CYP2C9", "DPYD"}
        if total <= 0.5:
            return PM
        return IM if total <= 1.5 else NM

    @pytest.mark.parametrize("gene", ["CYP2C19", "SLCO1B1", "TPMT", "NUDT15", "CYP3A5"])
    def test_function_defined_genes(self, gene: str) -> None:
        classes = self._classes(gene)
        for first, second in itertools.combinations_with_replacement(sorted(classes), 2):
            expected = self._by_function(gene, classes[first], classes[second])
            actual = classify_phenotype(f"{first}/{second}", gene)["phenotype"]
            assert actual == expected, f"{gene} {first}/{second}"

    @pytest.mark.parametrize("gene", ["CYP2D6", "CYP2C9", "DPYD"])
    def test_activity_score_defined_genes(self, gene: str) -> None:
        values = {name: d["activity_value"] for name, d in _BUILTIN_ALLELE_DEFINITIONS[gene].items()}
        for first, second in itertools.combinations_with_replacement(sorted(values), 2):
            expected = self._by_activity_score(gene, values[first] + values[second])
            actual = classify_phenotype(f"{first}/{second}", gene)["phenotype"]
            assert actual == expected, f"{gene} {first}/{second}"


class TestSingleSourceOfTruth:
    """The metabolism front end holds no table of its own."""

    def test_allele_function_table_is_the_canonical_table(self) -> None:
        assert default_allele_function_table() == {gene: dict(t) for gene, t in _ACTIVITY_SCORE_TABLES.items()}

    def test_allele_function_table_is_an_independent_copy(self) -> None:
        table = default_allele_function_table()
        table["CYP2D6"]["*4"] = 99.0
        assert _ACTIVITY_SCORE_TABLES["CYP2D6"]["*4"] == 0.0
        assert default_allele_function_table()["CYP2D6"]["*4"] == 0.0

    def test_cyp2d6_star10_uses_the_2019_consensus_value(self) -> None:
        # The metabolism table used the pre-2019 value 0.5.
        assert default_allele_function_table()["CYP2D6"]["*10"] == 0.25
        assert predict_metabolizer({"diplotype": "*4/*10"}, "CYP2D6")["activity_score"] == 0.25

    def test_every_calling_definition_carries_the_canonical_activity_value(self) -> None:
        for gene, alleles in _BUILTIN_ALLELE_DEFINITIONS.items():
            for name, definition in alleles.items():
                assert definition["activity_value"] == _ACTIVITY_SCORE_TABLES[gene][name], f"{gene} {name}"

    def test_both_apis_agree_on_every_diplotype_of_every_gene(self) -> None:
        checked = 0
        for gene, table in _ACTIVITY_SCORE_TABLES.items():
            for first, second in itertools.combinations_with_replacement(sorted(table), 2):
                diplotype = f"{first}/{second}"
                classified = classify_phenotype(diplotype, gene)
                predicted = predict_metabolizer({"diplotype": diplotype}, gene)
                assert predicted["activity_score"] == classified["activity_score"], f"{gene} {diplotype}"
                assert predicted["phenotype"] == _token(classified["phenotype"]), f"{gene} {diplotype}"
                checked += 1
        assert checked > 400

    def test_classify_metabolizer_is_the_phenotype_module_mapping(self) -> None:
        for gene in _PHENOTYPE_THRESHOLDS:
            for score in (0.0, 0.25, 0.5, 1.0, 1.25, 1.5, 2.0, 2.25, 2.5, 3.0):
                assert classify_metabolizer(score, gene) == _token(predict_metabolizer_status(score, gene))

    def test_thresholds_reported_by_the_front_end_are_the_ones_applied(self) -> None:
        reported = predict_metabolizer({"diplotype": "*1/*4"}, "CYP2D6")["thresholds_used"]
        assert reported["intermediate"]["high"] == 1.25
        assert reported["normal"]["low"] == 1.25 and reported["normal"]["low_inclusive"] is True
        assert reported["ultrarapid"]["low"] == 2.25 and reported["ultrarapid"]["low_inclusive"] is False

    def test_front_end_covers_every_gene_with_alleles(self) -> None:
        # SLCO1B1 and NUDT15 had no table in the metabolism module.
        assert predict_metabolizer({"diplotype": "*1/*5"}, "SLCO1B1")["phenotype"] == "intermediate"
        assert predict_metabolizer({"diplotype": "*2/*3"}, "NUDT15")["phenotype"] == "poor"


class TestScoreBands:
    """The score -> phenotype mapping covers every score, without gaps or overlaps."""

    @pytest.mark.parametrize("gene", sorted(_PHENOTYPE_THRESHOLDS))
    def test_bands_tile_from_zero_to_infinity(self, gene: str) -> None:
        bands = get_phenotype_thresholds(gene)
        assert bands[0]["lower_bound"] == 0.0 and bands[0]["lower_inclusive"] is True
        assert bands[-1]["upper_bound"] == float("inf") and bands[-1]["upper_inclusive"] is True
        for previous, current in itertools.pairwise(bands):
            assert current["lower_bound"] == previous["upper_bound"]
            # exactly one of the two bands owns the shared boundary score
            assert current["lower_inclusive"] != previous["upper_inclusive"]
        assert all(band["gene_specific"] for band in bands)

    @pytest.mark.parametrize("gene", [*sorted(_PHENOTYPE_THRESHOLDS), "NOTAGENE"])
    def test_phenotype_never_decreases_as_activity_rises(self, gene: str) -> None:
        order = [PM, IM, NM, RM, UM]
        ranks = [order.index(predict_metabolizer_status(i * 0.05, gene)) for i in range(0, 121)]
        assert ranks == sorted(ranks), f"{gene}: a score between two bands fell into the wrong one"

    def test_cyp2d6_boundaries(self) -> None:
        expected = {
            0.0: PM,
            0.1: IM,  # between 0 and 0.25: used to fall through to ultrarapid
            0.25: IM,
            1.0: IM,
            1.1: IM,
            1.25: NM,
            2.0: NM,
            2.25: NM,  # upper end of normal is inclusive
            2.3: UM,
            2.5: UM,
            3.0: UM,
        }
        for score, phenotype in expected.items():
            assert predict_metabolizer_status(score, "CYP2D6") == phenotype, score

    def test_cyp2c19_scores_use_the_allele_scale(self) -> None:
        # *1 counts 1.0 and *17 counts 1.5, so *1/*1 is 2.0 and *1/*17 is 2.5.
        expected = {0.0: PM, 1.0: IM, 1.5: IM, 2.0: NM, 2.5: RM, 3.0: UM}
        for score, phenotype in expected.items():
            assert predict_metabolizer_status(score, "CYP2C19") == phenotype, score
            assert classify_metabolizer(score, "CYP2C19") == _token(phenotype), score

    def test_slco1b1_scores(self) -> None:
        # *5 and *15 are no function (0): *5/*5 = 0 poor, *1/*5 = 1.0 decreased, *1/*1 = 2.0 normal.
        expected = {0.0: PM, 1.0: IM, 2.0: NM}
        for score, phenotype in expected.items():
            assert predict_metabolizer_status(score, "SLCO1B1") == phenotype, score

    @pytest.mark.parametrize("bad_score", [-0.25, -1.0, float("nan")])
    def test_invalid_scores_are_rejected_not_called_ultrarapid(self, bad_score: float) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            predict_metabolizer_status(bad_score, "CYP2D6")
        with pytest.raises(ValueError, match="non-negative"):
            classify_metabolizer(bad_score, "CYP2D6")

    def test_unknown_gene_reports_the_default_bands_it_applies(self) -> None:
        bands = get_phenotype_thresholds("NOTAGENE")
        reference = get_phenotype_thresholds("CYP2D6")
        assert [(b["lower_bound"], b["upper_bound"], b["phenotype"]) for b in bands] == [
            (b["lower_bound"], b["upper_bound"], b["phenotype"]) for b in reference
        ]
        assert all(b["gene_specific"] is False and "NOTAGENE" in b["note"] for b in bands)
        assert predict_metabolizer_status(1.0, "NOTAGENE") == IM


class TestFunctionRulesVersusSums:
    """CYP2C19 is defined by allele function, which a sum of activity values cannot always express."""

    def test_decreased_function_pair_that_sums_like_star1_star1_is_not_normal(self) -> None:
        assert determine_diplotype("*9", "*17", "CYP2C19").activity_score == 2.0
        assert determine_diplotype("*1", "*1", "CYP2C19").activity_score == 2.0
        assert classify_phenotype("*9/*17", "CYP2C19")["phenotype"] == IM
        assert classify_phenotype("*1/*1", "CYP2C19")["phenotype"] == NM

    def test_score_only_mapping_cannot_tell_those_pairs_apart(self) -> None:
        # Documented limitation of predict_metabolizer_status(score, gene): 2.0 reads as normal.
        assert predict_metabolizer_status(2.0, "CYP2C19") == NM

    def test_diplotype_scored_with_a_custom_table_is_classified_from_its_stored_score(self) -> None:
        custom = {"*1": 1.0, "*2": 1.0}  # a table in which *2 is a normal-function allele
        diplotype = determine_diplotype("*1", "*2", "CYP2C19", scoring_table=custom)
        assert diplotype.activity_score == 2.0
        assert classify_phenotype(diplotype, "CYP2C19")["phenotype"] == NM

    def test_increased_function_slco1b1_allele_has_no_phenotype_and_raises(self) -> None:
        table = {"*1": 1.0, "*X": 1.5}
        with pytest.raises(ValueError, match="No SLCO1B1 phenotype"):
            predict_metabolizer({"diplotype": "*1/*X"}, "SLCO1B1", activity_scores=table)

    def test_ugt1a1_star37_with_star1_is_intermediate_whether_recorded_as_decreased_or_no_function(self) -> None:
        # CPIC (Gammal 2016) calls *1/*37 intermediate. Going by function class keeps that true
        # for either recorded value of *37 (0.5 decreased, 0.0 no function).
        for value in (0.0, 0.5):
            table = {"*1": 1.0, "*37": value}
            result = predict_metabolizer({"diplotype": "*1/*37"}, "UGT1A1", activity_scores=table)
            assert result["phenotype"] == "intermediate", value


class TestPopulationFrequencies:
    def test_cyp2c19_hardy_weinberg_distribution(self) -> None:
        freqs = population_phenotype_frequencies({"*1": 0.5, "*17": 0.5}, "CYP2C19")
        assert freqs == {"NM": 0.25, "RM": 0.5, "UM": 0.25}

    def test_cyp2d6_hardy_weinberg_distribution(self) -> None:
        freqs = population_phenotype_frequencies({"*1": 0.7, "*4": 0.3}, "CYP2D6")
        assert freqs["NM"] == pytest.approx(0.49)
        assert freqs["IM"] == pytest.approx(0.42)  # *1/*4 has activity score 1.0
        assert freqs["PM"] == pytest.approx(0.09)

    def test_slco1b1_distribution(self) -> None:
        freqs = population_phenotype_frequencies({"*1": 0.8, "*5": 0.2}, "SLCO1B1")
        assert freqs == {"NM": 0.64, "IM": 0.32, "PM": 0.04}


class TestPhenotypeLabels:
    @pytest.mark.parametrize(
        ("label", "expected"),
        [
            ("Poor Metabolizer", PM),
            ("poor", PM),
            ("PM", PM),
            ("  Poor  ", PM),
            ("Intermediate Metabolizer", IM),
            ("im", IM),
            ("Normal Metabolizer", NM),
            ("Extensive Metabolizer", NM),
            ("EM", NM),
            ("Rapid Metabolizer", RM),
            ("rapid", RM),
            ("Ultrarapid Metabolizer", UM),
            ("ultra-rapid", UM),
            ("Ultra Rapid Metabolizer", UM),
            ("UM", UM),
            ("Indeterminate", MetabolizerPhenotype.INDETERMINATE),
            (MetabolizerPhenotype.INTERMEDIATE, IM),
        ],
    )
    def test_from_string(self, label: str | MetabolizerPhenotype, expected: MetabolizerPhenotype) -> None:
        assert MetabolizerPhenotype.from_string(label) is expected

    @pytest.mark.parametrize("label", ["superrapid", "", "poor-ish", "metabolizer", "None"])
    def test_from_string_rejects_unrecognised_labels(self, label: str) -> None:
        with pytest.raises(ValueError, match="Unrecognised metabolizer phenotype"):
            MetabolizerPhenotype.from_string(label)


class TestDoseAdjustmentPhenotypeLabels:
    POOR_FORMS = ["poor", "Poor", " POOR ", "Poor Metabolizer", "poor metabolizer", "PM", "pm", PM]

    @pytest.mark.parametrize("label", POOR_FORMS, ids=[str(f) for f in POOR_FORMS])
    def test_poor_metabolizer_in_any_form_finds_the_poor_guideline(self, label: str | MetabolizerPhenotype) -> None:
        result = dose_adjustment(label, "codeine")
        assert result == dose_adjustment("poor", "codeine")
        assert result["adjusted_dose_fraction"] == 0.0
        assert result["evidence_level"] == "A"
        assert result["guideline_source"] == "CPIC"
        assert result["metabolizer_status"] == "poor"

    @pytest.mark.parametrize("label", ["Ultrarapid Metabolizer", "ultra-rapid", "UM", "ultrarapid"])
    def test_ultrarapid_in_any_form(self, label: str) -> None:
        result = dose_adjustment(label, "codeine")
        assert result["metabolizer_status"] == "ultrarapid"
        assert "morphine" in result["recommendation"]

    def test_full_name_finds_a_dose_reduction(self) -> None:
        result = dose_adjustment("Intermediate Metabolizer", "azathioprine")
        assert result["adjusted_dose_fraction"] == 0.5

    def test_custom_guidelines_may_key_phenotypes_by_full_name_or_abbreviation(self) -> None:
        custom = {
            "exampledrug": {
                "Poor Metabolizer": {
                    "recommendation": "Halve the dose",
                    "adjusted_dose_fraction": 0.5,
                    "evidence_level": "B",
                    "guideline_source": "local",
                },
                "NM": {
                    "recommendation": "Standard dose",
                    "adjusted_dose_fraction": 1.0,
                    "evidence_level": "B",
                    "guideline_source": "local",
                },
            }
        }
        poor = dose_adjustment("poor", "ExampleDrug", dose_guidelines=custom)
        assert poor["recommendation"] == "Halve the dose"
        normal = dose_adjustment("Normal Metabolizer", "exampledrug", dose_guidelines=custom)
        assert normal["adjusted_dose_fraction"] == 1.0

    def test_unrecognised_status_still_falls_back_to_no_guideline(self) -> None:
        result = dose_adjustment("superrapid", "codeine")
        assert result["evidence_level"] == "D"
        assert result["metabolizer_status"] == "superrapid"
        assert "superrapid" in result["recommendation"]
        assert dose_adjustment(MetabolizerPhenotype.INDETERMINATE, "codeine")["evidence_level"] == "D"
