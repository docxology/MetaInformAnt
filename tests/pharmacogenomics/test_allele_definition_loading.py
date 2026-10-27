"""Allele definition loading: the ``source`` / ``filepath`` contract and genes without definitions.

``load_allele_definitions`` used to ignore ``filepath`` unless ``source="file"`` and to fall back
to the built-in table whenever ``filepath`` was None, whatever ``source`` said, so a custom
file could be silently dropped. UGT1A1 has activity scores but no calling definitions; that gap
must be stated, not look like a typo.

Definition files are real JSON/TSV files written to ``tmp_path``.

REAL IMPLEMENTATION -- all tests use real implementations.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.pharmacogenomics.alleles import star_allele
from metainformant.pharmacogenomics.alleles.diplotype import (
    _ACTIVITY_SCORE_TABLES,
    determine_diplotype,
    phased_diplotype,
)
from metainformant.pharmacogenomics.alleles.phenotype import MetabolizerPhenotype, classify_phenotype
from metainformant.pharmacogenomics.alleles.star_allele import (
    _BUILTIN_ALLELE_DEFINITIONS,
    call_star_alleles,
    load_allele_definitions,
)


def _definition(variants: list[str], function: str, activity: float) -> dict[str, object]:
    return {"defining_variants": variants, "function": function, "activity_value": activity}


@pytest.fixture()
def custom_json(tmp_path: Path) -> Path:
    """A single-gene definition file whose alleles are not in any built-in table."""
    path = tmp_path / "custom_defs.json"
    path.write_text(
        json.dumps(
            {
                "*1": _definition([], "Normal function", 1.0),
                "*77": _definition(["rs7700001"], "No function", 0.0),
                "*78": _definition(["rs7800001", "rs7800002"], "Decreased function", 0.5),
            }
        ),
        encoding="utf-8",
    )
    return path


CUSTOM_NAMES = ["*1", "*77", "*78"]


class TestSourceAndFilepathContract:
    def test_filepath_alone_loads_the_file(self, custom_json: Path) -> None:
        alleles = load_allele_definitions("CYP2E1", filepath=custom_json)
        assert [a.name for a in alleles] == CUSTOM_NAMES
        assert all(a.gene == "CYP2E1" for a in alleles)

    def test_filepath_alone_replaces_a_built_in_gene_instead_of_being_ignored(self, custom_json: Path) -> None:
        alleles = load_allele_definitions("CYP2D6", filepath=custom_json)
        assert [a.name for a in alleles] == CUSTOM_NAMES  # not the 12 built-in CYP2D6 alleles

    def test_source_file_with_filepath_loads_the_file(self, custom_json: Path) -> None:
        alleles = load_allele_definitions("CYP2E1", source="file", filepath=custom_json)
        assert [a.name for a in alleles] == CUSTOM_NAMES
        by_name = {a.name: a for a in alleles}
        assert by_name["*78"].defining_variants == frozenset({"rs7800001", "rs7800002"})
        assert by_name["*78"].activity_value == 0.5
        assert by_name["*77"].is_no_function

    def test_gene_keyed_json_file(self, tmp_path: Path) -> None:
        path = tmp_path / "multi_gene.json"
        path.write_text(
            json.dumps(
                {
                    "CYP2E1": {"*1": _definition([], "Normal function", 1.0)},
                    "CYP2B6": {"*6": _definition(["rs3745274"], "Decreased function", 0.5)},
                }
            ),
            encoding="utf-8",
        )
        assert [a.name for a in load_allele_definitions("cyp2b6", filepath=path)] == ["*6"]

    def test_tsv_file_with_filepath_alone(self, tmp_path: Path) -> None:
        path = tmp_path / "defs.tsv"
        path.write_text(
            "allele\tdefining_variants\tfunction\tactivity_value\n"
            "*1\t\tNormal function\t1.0\n"
            "*5\trs5000001;rs5000002\tNo function\t0.0\n",
            encoding="utf-8",
        )
        alleles = load_allele_definitions("CYP2E1", filepath=path)
        assert [a.name for a in alleles] == ["*1", "*5"]
        assert alleles[1].defining_variants == frozenset({"rs5000001", "rs5000002"})

    def test_missing_file_raises_file_not_found(self, tmp_path: Path) -> None:
        missing = tmp_path / "absent.json"
        with pytest.raises(FileNotFoundError):
            load_allele_definitions("CYP2D6", filepath=missing)
        with pytest.raises(FileNotFoundError):
            load_allele_definitions("CYP2D6", source="file", filepath=missing)

    def test_source_file_without_filepath_raises_instead_of_loading_the_built_in_table(self) -> None:
        with pytest.raises(ValueError, match="requires a filepath"):
            load_allele_definitions("CYP2D6", source="file")

    def test_source_cpic_with_filepath_raises_instead_of_ignoring_the_file(self, custom_json: Path) -> None:
        with pytest.raises(ValueError, match="only used with source='file'"):
            load_allele_definitions("CYP2D6", source="cpic", filepath=custom_json)

    @pytest.mark.parametrize("source", ["bogus", "builtin", "files", ""])
    def test_invalid_source_raises_with_or_without_a_filepath(self, source: str, custom_json: Path) -> None:
        with pytest.raises(ValueError, match="Invalid source"):
            load_allele_definitions("CYP2D6", source=source)
        with pytest.raises(ValueError, match="Invalid source"):
            load_allele_definitions("CYP2D6", source=source, filepath=custom_json)

    def test_source_is_case_insensitive(self, custom_json: Path) -> None:
        assert [a.name for a in load_allele_definitions("CYP2E1", source="File", filepath=custom_json)] == CUSTOM_NAMES
        built_in = [a.name for a in load_allele_definitions("CYP2D6")]
        assert [a.name for a in load_allele_definitions("CYP2D6", source=" CPIC ")] == built_in

    def test_default_and_explicit_cpic_load_the_built_in_table(self) -> None:
        built_in = [a.name for a in load_allele_definitions("CYP2D6")]
        assert built_in == list(_BUILTIN_ALLELE_DEFINITIONS["CYP2D6"])
        assert [a.name for a in load_allele_definitions("CYP2D6", source="cpic")] == built_in
        assert [a.name for a in load_allele_definitions("CYP2D6", source=None, filepath=None)] == built_in

    def test_loaded_file_definitions_drive_allele_calling(self, custom_json: Path) -> None:
        definitions = load_allele_definitions("CYP2E1", filepath=custom_json)
        called = call_star_alleles({"rs7800001", "rs7800002"}, "CYP2E1", allele_definitions=definitions)
        assert [a.name for a in called] == ["*78"]


class TestGenesWithoutCallingDefinitions:
    def test_ugt1a1_error_states_the_gap(self) -> None:
        with pytest.raises(ValueError) as caught:
            load_allele_definitions("UGT1A1")
        message = str(caught.value)
        assert "UGT1A1" in message
        assert "no built-in star-allele definitions" in message
        assert "TA-repeat" in message
        assert "determine_diplotype()" in message and "source='file'" in message
        assert "not found" not in message  # it is a known gene, not a typo

    def test_unknown_gene_error_is_unchanged(self) -> None:
        with pytest.raises(ValueError, match="not found in built-in definitions") as caught:
            load_allele_definitions("NOTAGENE")
        assert "Available genes" in str(caught.value)

    def test_calling_and_phasing_report_the_same_gap(self) -> None:
        with pytest.raises(ValueError, match="no built-in star-allele definitions"):
            call_star_alleles({"rs4148323"}, "UGT1A1")
        with pytest.raises(ValueError, match="no built-in star-allele definitions"):
            phased_diplotype({"rs4148323"}, {"rs4148323": 0}, "UGT1A1")

    def test_registry_lists_every_gene_that_has_scores_but_no_definitions(self) -> None:
        registry = star_allele._GENES_WITHOUT_BUILTIN_DEFINITIONS
        without_definitions = set(_ACTIVITY_SCORE_TABLES) - set(_BUILTIN_ALLELE_DEFINITIONS)
        assert without_definitions == set(registry)
        assert "UGT1A1" in without_definitions
        assert all(registry[gene].strip() for gene in without_definitions)

    def test_registered_genes_have_no_built_in_definitions(self) -> None:
        registry = star_allele._GENES_WITHOUT_BUILTIN_DEFINITIONS
        assert not set(registry) & set(_BUILTIN_ALLELE_DEFINITIONS)

    def test_every_defined_allele_has_an_activity_score(self) -> None:
        for gene, alleles in _BUILTIN_ALLELE_DEFINITIONS.items():
            assert set(alleles) <= set(_ACTIVITY_SCORE_TABLES[gene]), gene

    def test_ugt1a1_diplotypes_can_still_be_scored_and_classified_from_allele_names(self) -> None:
        assert determine_diplotype("*1", "*28", "UGT1A1").activity_score == 1.5
        assert classify_phenotype("*1/*28", "UGT1A1")["phenotype"] == MetabolizerPhenotype.INTERMEDIATE
        assert classify_phenotype("*28/*28", "UGT1A1")["phenotype"] == MetabolizerPhenotype.POOR

    def test_ugt1a1_calling_works_through_definitions_supplied_from_a_file(self, tmp_path: Path) -> None:
        # The supported route named in the error message. The TA-repeat tokens are this test's
        # own convention (repeat length is not a presence/absence variant); they are not shipped.
        path = tmp_path / "ugt1a1.json"
        path.write_text(
            json.dumps(
                {
                    "*1": _definition([], "Normal function", 1.0),
                    "*6": _definition(["rs4148323"], "Decreased function", 0.5),
                    "*28": _definition(["UGT1A1_TA7"], "Decreased function", 0.5),
                }
            ),
            encoding="utf-8",
        )
        definitions = load_allele_definitions("UGT1A1", filepath=path)
        called = call_star_alleles({"UGT1A1_TA7"}, "UGT1A1", allele_definitions=definitions)
        assert [a.name for a in called] == ["*28"]
        diplotype = determine_diplotype(called[0], "*28", "UGT1A1")
        assert classify_phenotype(diplotype, "UGT1A1")["phenotype"] == MetabolizerPhenotype.POOR
