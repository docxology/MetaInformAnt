"""Tests for AntWiki record handling: AntWikiRecord, filtering, matrices, reports.

Complements test_phenotype_basic.py / test_phenotype_comprehensive.py which
cover load_antwiki_json; here the AntWikiRecord object API is exercised with
real deterministic data.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.phenotype.data.antwiki import (
    AntWikiRecord,
    create_phenotype_matrix,
    filter_antwiki_records,
    find_similar_species,
    generate_antwiki_report,
    get_phenotype_distribution,
    save_antwiki_json,
)


def _record(species: str, genus: str, **kwargs) -> AntWikiRecord:
    data = {"species_name": species, "genus": genus}
    data.update(kwargs)
    return AntWikiRecord(data)


@pytest.fixture
def records() -> list[AntWikiRecord]:
    return [
        _record(
            "pennsylvanicus",
            "Camponotus",
            subfamily="Formicinae",
            tribe="Camponotini",
            morphology={"body_length_mm": 12.0, "head_width_mm": 2.0},
            behavior={"foraging_strategy": "generalist"},
            confidence_score=0.9,
        ),
        _record(
            "rufa",
            "Formica",
            subfamily="Formicinae",
            morphology={"body_length_mm": 8.0, "head_width_mm": 1.5},
            behavior={"foraging_strategy": "generalist"},
            confidence_score=0.7,
        ),
        _record(
            "niger",
            "Lasius",
            subfamily="Formicinae",
            morphology={"body_length_mm": 3.0},
            confidence_score=0.4,
        ),
    ]


class TestAntWikiRecord:
    def test_basic_fields_and_taxonomic_path(self, records):
        rec = records[0]
        assert rec.species_name == "pennsylvanicus"
        assert rec.genus == "Camponotus"
        assert rec.subfamily == "Formicinae"
        assert rec.tribe == "Camponotini"
        assert rec.taxonomic_path == ["Camponotini", "Formicinae", "Camponotus", "pennsylvanicus"]

    def test_taxonomic_path_skips_missing_ranks(self, records):
        rec = _record("barbatus", "Pogonomyrmex")
        assert rec.taxonomic_path == ["Pogonomyrmex", "barbatus"]

    def test_missing_species_name_raises(self):
        with pytest.raises(ValueError, match="species_name"):
            AntWikiRecord({"genus": "Camponotus"})

    def test_missing_genus_raises(self):
        with pytest.raises(ValueError, match="genus"):
            AntWikiRecord({"species_name": "pennsylvanicus"})

    def test_phenotype_extraction(self, records):
        rec = records[0]
        assert rec.phenotypes == {
            "body_length": 12.0,
            "head_width": 2.0,
            "foraging_strategy": "generalist",
        }
        assert rec.has_phenotype("body_length")
        assert not rec.has_phenotype("nest_type")
        assert rec.get_phenotype_value("head_width") == 2.0
        assert rec.get_phenotype_value("nonexistent") is None

    def test_raw_data_and_to_dict_roundtrip(self, records):
        rec = records[0]
        d = rec.to_dict()
        assert d["species_name"] == "pennsylvanicus"
        assert d["confidence_score"] == 0.9
        rebuilt = AntWikiRecord(
            {
                "species_name": d["species_name"],
                "genus": d["genus"],
                "subfamily": d["subfamily"],
                "tribe": d["tribe"],
                "morphology": d["morphology"],
                "behavior": d["behavior"],
                "confidence_score": d["confidence_score"],
            }
        )
        assert rebuilt.phenotypes == rec.phenotypes


class TestFilterRecords:
    def test_filter_by_genus(self, records):
        out = filter_antwiki_records(records, genus="Formica")
        assert [r.species_name for r in out] == ["rufa"]

    def test_filter_by_subfamily(self, records):
        out = filter_antwiki_records(records, subfamily="Formicinae")
        assert len(out) == 3

    def test_filter_by_min_confidence(self, records):
        out = filter_antwiki_records(records, min_confidence=0.5)
        assert [r.species_name for r in out] == ["pennsylvanicus", "rufa"]

    def test_filter_by_required_phenotypes(self, records):
        out = filter_antwiki_records(records, required_phenotypes=["head_width", "foraging_strategy"])
        assert [r.species_name for r in out] == ["pennsylvanicus", "rufa"]

    def test_no_filters_returns_all(self, records):
        assert filter_antwiki_records(records) == records


class TestPhenotypeDistribution:
    def test_numeric_distribution(self, records):
        stats = get_phenotype_distribution(records, "body_length")
        assert stats["values_found"] == 3
        assert stats["coverage"] == 1.0
        assert stats["unique_values"] == 3
        assert stats["mean"] == pytest.approx((12.0 + 8.0 + 3.0) / 3)
        assert stats["min"] == 3.0
        assert stats["max"] == 12.0
        assert stats["median"] == 8.0

    def test_absent_phenotype(self, records):
        stats = get_phenotype_distribution(records, "nonexistent")
        assert stats["values_found"] == 0
        assert "mean" not in stats


class TestFindSimilarSpecies:
    def test_matching_foraging_and_morphology_ranks_first(self, records):
        target = records[0]
        out = find_similar_species(records, target, top_k=2)
        assert all(species != target.species_name for species, _ in out)
        # Similarity scores sorted descending
        scores = [s for _, s in out]
        assert scores == sorted(scores, reverse=True)
        # rufa shares foraging_strategy and is closer in size than niger
        assert out[0][0].species_name == "rufa"

    def test_top_k_limits_results(self, records):
        out = find_similar_species(records, records[0], top_k=1)
        assert len(out) == 1


class TestPhenotypeMatrix:
    def test_matrix_shape_and_values(self, records):
        species, phenotypes, matrix = create_phenotype_matrix(records, ["body_length"])
        assert species == ["Camponotus pennsylvanicus", "Formica rufa", "Lasius niger"]
        assert phenotypes == ["body_length"]
        assert [row[0] for row in matrix] == [12.0, 8.0, 3.0]

    def test_all_phenotypes_and_missing_cells(self, records):
        species, phenotypes, matrix = create_phenotype_matrix(records)
        assert "body_length" in phenotypes
        assert "foraging_strategy" in phenotypes
        i_bl = phenotypes.index("body_length")
        i_fw = phenotypes.index("foraging_strategy")
        # niger has no foraging strategy
        assert matrix[2][i_fw] is None
        assert matrix[0][i_bl] == 12.0

    def test_empty_records(self):
        assert create_phenotype_matrix([]) == ([], [], [])


class TestSaveAndReport:
    def test_save_antwiki_json(self, records, tmp_path: Path):
        out = tmp_path / "records.json"
        save_antwiki_json(records, out)
        payload = json.loads(out.read_text())
        assert payload["metadata"]["total_records"] == 3
        assert payload["metadata"]["data_source"] == "antwiki"
        assert len(payload["records"]) == 3
        assert payload["records"][0]["species_name"] == "pennsylvanicus"

    def test_report_contents(self, records, tmp_path: Path):
        out = tmp_path / "report.txt"
        report = generate_antwiki_report(records, output_path=out)
        assert out.exists()
        assert "ANTWIKI DATA SUMMARY REPORT" in report
        assert "Total Records: 3" in report
        assert "Camponotus: 1 species" in report
        assert "body_length: 3/3 (100.0%)" in report
        assert "Average Confidence Score: 0.67" in report

    def test_report_without_output_path(self, records):
        report = generate_antwiki_report(records)
        assert "Total Records: 3" in report
