"""Zero-mocks tests for gwas.data.traits and gwas.data.extraction.

Real CSV files and real ZIP archives written to tmp_path; no fake loaders.
"""

from __future__ import annotations

import zipfile

import pytest

from metainformant.gwas.data import extraction, traits


class TestLoadTraits:
    def test_loads_last_non_id_column(self, tmp_path) -> None:
        f = tmp_path / "pheno.csv"
        f.write_text("sample_id,body_weight\nS1,10.5\nS2,11.0\nS3,notanumber\nS4,9.25\n")
        vals = traits.load_traits(f)
        assert vals == pytest.approx([10.5, 11.0, 9.25])

    def test_explicit_trait_column(self, tmp_path) -> None:
        f = tmp_path / "pheno.csv"
        f.write_text("sample_id,aggression,body_weight\nS1,1.5,10.5\nS2,2.5,11.0\n")
        vals = traits.load_traits(f, trait_column="aggression")
        assert vals == pytest.approx([1.5, 2.5])

    def test_empty_file(self, tmp_path) -> None:
        f = tmp_path / "empty.csv"
        f.write_text("")
        assert traits.load_traits(f) == []

    def test_header_only(self, tmp_path) -> None:
        f = tmp_path / "hdr.csv"
        f.write_text("sample_id,trait\n")
        assert traits.load_traits(f) == []


class TestLoadTraitsWithIds:
    def test_pairs_ids_and_values(self, tmp_path) -> None:
        f = tmp_path / "pheno.csv"
        f.write_text("sample_id,trait_value\nS1,1.0\nS2,2.0\nS3,bad\nS4,4.0\n")
        ids, vals = traits.load_traits_with_ids(f)
        assert ids == ["S1", "S2", "S4"]
        assert vals == pytest.approx([1.0, 2.0, 4.0])

    def test_auto_detect_id_column(self, tmp_path) -> None:
        f = tmp_path / "pheno.csv"
        f.write_text("colony,wing_mm\nC15,0.9\nC16,0.95\n")
        ids, vals = traits.load_traits_with_ids(f, id_column="colony", trait_column="wing_mm")
        assert ids == ["C15", "C16"]
        assert vals == pytest.approx([0.9, 0.95])

    def test_missing_id_column_returns_empty(self, tmp_path) -> None:
        f = tmp_path / "pheno.csv"
        f.write_text("a,b\n1,2\n")
        ids, vals = traits.load_traits_with_ids(f, id_column="nonexistent")
        assert ids == [] and vals == []


class TestParseFastqFilename:
    @pytest.mark.parametrize(
        "name,expected",
        [
            ("C15G_R1.fastq", {"colony": "C15", "biological_group_code": "G", "read": "R1", "sample_id": "C15G"}),
            (
                "M03WORK_R2.fastq",
                {"colony": "M03", "biological_group_code": "WORK", "read": "R2", "sample_id": "M03WORK"},
            ),
            (
                "R21ITQ_R1.fastq-001.gz",
                {"colony": "R21", "biological_group_code": "ITQ", "read": "R1", "sample_id": "R21ITQ"},
            ),
            (
                "I44IV_R2.fastq-12.gz",
                {"colony": "I44", "biological_group_code": "IV", "read": "R2", "sample_id": "I44IV"},
            ),
        ],
    )
    def test_valid_names(self, name: str, expected: dict) -> None:
        assert extraction.parse_fastq_filename(name) == expected

    @pytest.mark.parametrize("name", ["notes.txt", "X15Q_R1.fastq", "C15_R1.fastq", "C15G_R3.fastq", ""])
    def test_invalid_names(self, name: str) -> None:
        assert extraction.parse_fastq_filename(name) is None


class TestLineageFromColony:
    def test_known_lineages(self) -> None:
        assert extraction.lineage_from_colony("C15") == "unknown"  # C not in the I/M/R mapping
        assert extraction.lineage_from_colony("I12") == "I"
        assert extraction.lineage_from_colony("M03") == "M"
        assert extraction.lineage_from_colony("R21") == "R"

    def test_unknown_prefix(self) -> None:
        assert extraction.lineage_from_colony("Z99") == "unknown"


class TestBuildSampleFileMap:
    def test_maps_gz_files(self, tmp_path) -> None:
        (tmp_path / "C15G_R1.fastq-001.gz").write_bytes(b"x")
        (tmp_path / "C15G_R2.fastq-001.gz").write_bytes(b"x")
        (tmp_path / "I12ITQ_R1.fastq-001.gz").write_bytes(b"x")
        (tmp_path / "unrelated.txt").write_text("skip me")
        smap = extraction.build_sample_file_map(tmp_path)
        assert set(smap.keys()) == {"C15G", "I12ITQ"}
        assert smap["C15G"]["R1"]["source_type"] == "standalone"
        assert smap["C15G"]["R1"]["colony"] == "C15"

    def test_maps_zip_contents(self, tmp_path) -> None:
        zpath = tmp_path / "batch1.zip"
        with zipfile.ZipFile(zpath, "w") as zf:
            zf.writestr("C15G_R1.fastq", "@read\nACGT\n")
            zf.writestr("C15G_R2.fastq", "@read\nACGT\n")
        smap = extraction.build_sample_file_map(tmp_path)
        assert "C15G" in smap
        assert smap["C15G"]["R1"]["source_type"] == "zip"
        assert smap["C15G"]["R1"]["inner_filename"] == "C15G_R1.fastq"
        assert smap["C15G"]["R2"]["source_type"] == "zip"

    def test_accessions_filter(self, tmp_path) -> None:
        (tmp_path / "C15G_R1.fastq-001.gz").write_bytes(b"x")
        (tmp_path / "I12ITQ_R1.fastq-001.gz").write_bytes(b"x")
        smap = extraction.build_sample_file_map(tmp_path, accessions={"C15G"})
        assert set(smap.keys()) == {"C15G"}

    def test_missing_dir(self, tmp_path) -> None:
        assert extraction.build_sample_file_map(tmp_path / "nope") == {}

    def test_bad_zip_reported_not_raised(self, tmp_path) -> None:
        bad = tmp_path / "bad.zip"
        bad.write_bytes(b"not a zip at all")
        (tmp_path / "C15G_R1.fastq-001.gz").write_bytes(b"x")
        smap = extraction.build_sample_file_map(tmp_path)
        assert set(smap.keys()) == {"C15G"}  # bad ZIP skipped, good data survives
