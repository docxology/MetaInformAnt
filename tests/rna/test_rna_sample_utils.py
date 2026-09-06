"""
Unit Tests for RNA sample utilities.

Verifies sample-ID extraction precedence, TSV metadata parsing, and
quantification-file discovery against real files on disk (tmp_path).
No mocks: every file-based behavior is exercised with real files.
"""

from pathlib import Path

from metainformant.rna.core.sample_utils import (
    extract_sample_id,
    find_quantification_file,
    quantification_file_candidates,
    read_sample_ids_from_metadata,
)


class TestExtractSampleId:
    def test_precedence_order(self) -> None:
        """Earlier SAMPLE_ID_COLUMNS entries win over later ones."""
        row = {
            "accession": "ACCESSION1",
            "sample_id": "SAMPLE1",
            "sra_run": "SRR2",
            "run": "SRR1",
        }
        assert extract_sample_id(row) == "SRR1"

    def test_case_insensitive_column_match(self) -> None:
        """Column names match case-insensitively against SAMPLE_ID_COLUMNS."""
        assert extract_sample_id({"Run": "SRR1"}) == "SRR1"
        assert extract_sample_id({"SRA_Run": "SRR2"}) == "SRR2"
        assert extract_sample_id({"Accession": "ACCESSION1"}) == "ACCESSION1"

    def test_skips_blank_and_whitespace_values(self) -> None:
        """Blank/whitespace-only values in a preferred column fall through."""
        row = {"run": "   ", "sample_id": None, "accession": "ACCESSION1"}
        assert extract_sample_id(row) == "ACCESSION1"

    def test_returns_empty_string_without_fallback(self) -> None:
        """No matching column and no fallback yields an empty string."""
        assert extract_sample_id({"tissue": "brain"}) == ""

    def test_returns_fallback_when_provided(self) -> None:
        """No matching column returns the caller-supplied fallback."""
        assert extract_sample_id({"tissue": "brain"}, fallback="FALLBACK1") == "FALLBACK1"

    def test_strips_surrounding_whitespace(self) -> None:
        """Matched values are returned stripped of surrounding whitespace."""
        assert extract_sample_id({"run": "  SRR1\t"}) == "SRR1"


class TestReadSampleIdsFromMetadata:
    def test_reads_ids_in_order_from_real_tsv(self, tmp_path: Path) -> None:
        """Sample IDs are parsed in row order from a real tab-delimited TSV."""
        metadata = tmp_path / "metadata.tsv"
        metadata.write_text(
            "run\ttissue\nSRR3\tbrain\nSRR1\twing\nSRR2\tantenna\n",
            encoding="utf-8",
        )
        assert read_sample_ids_from_metadata(metadata) == ["SRR3", "SRR1", "SRR2"]

    def test_header_only_returns_empty_list(self, tmp_path: Path) -> None:
        """A header-only TSV yields no sample IDs."""
        metadata = tmp_path / "metadata.tsv"
        metadata.write_text("run\ttissue\n", encoding="utf-8")
        assert read_sample_ids_from_metadata(metadata) == []


class TestQuantificationFileCandidates:
    def test_priority_order_with_sample_id(self, tmp_path: Path) -> None:
        """Candidates appear in priority order when a sample_id is given."""
        quant_dir = tmp_path / "SRR1"
        quant_dir.mkdir()
        candidates = quantification_file_candidates(quant_dir, sample_id="SRR1")
        assert candidates == [
            quant_dir / "abundance.tsv",
            quant_dir / "SRR1_abundance.tsv",
            quant_dir / "quant.sf",
        ]

    def test_no_per_sample_entry_without_sample_id(self, tmp_path: Path) -> None:
        """Without a sample_id the per-sample candidate is omitted."""
        quant_dir = tmp_path / "SRR1"
        quant_dir.mkdir()
        candidates = quantification_file_candidates(quant_dir)
        assert candidates == [
            quant_dir / "abundance.tsv",
            quant_dir / "quant.sf",
        ]


class TestFindQuantificationFile:
    def test_returns_first_existing_candidate(self, tmp_path: Path) -> None:
        """The highest-priority existing candidate is returned."""
        quant_dir = tmp_path / "SRR1"
        quant_dir.mkdir()
        (quant_dir / "abundance.tsv").write_text("target\tcount\n", encoding="utf-8")
        (quant_dir / "SRR1_abundance.tsv").write_text("later\n", encoding="utf-8")
        (quant_dir / "quant.sf").write_text("later_still\n", encoding="utf-8")
        assert find_quantification_file(quant_dir, sample_id="SRR1") == (
            quant_dir / "abundance.tsv"
        )

    def test_require_nonempty_skips_zero_byte_files(self, tmp_path: Path) -> None:
        """require_nonempty=True skips zero-byte candidates, returns non-empty."""
        quant_dir = tmp_path / "SRR1"
        quant_dir.mkdir()
        (quant_dir / "abundance.tsv").write_text("", encoding="utf-8")
        (quant_dir / "quant.sf").write_text("gene_id\tTPM\n", encoding="utf-8")
        assert find_quantification_file(quant_dir, sample_id="SRR1") == (
            quant_dir / "quant.sf"
        )

    def test_require_nonempty_false_returns_zero_byte_file(self, tmp_path: Path) -> None:
        """require_nonempty=False accepts a zero-byte candidate."""
        quant_dir = tmp_path / "SRR1"
        quant_dir.mkdir()
        empty = quant_dir / "abundance.tsv"
        empty.write_text("", encoding="utf-8")
        assert find_quantification_file(quant_dir, sample_id="SRR1", require_nonempty=False) == empty

    def test_glob_fallback_for_named_abundance_files(self, tmp_path: Path) -> None:
        """Named candidates missing: glob picks up <sample>_abundance.tsv."""
        quant_dir = tmp_path / "SRR1"
        quant_dir.mkdir()
        (quant_dir / "SRR1_abundance.tsv").write_text("target\tcount\n", encoding="utf-8")
        assert find_quantification_file(quant_dir, sample_id="SRR1") == (
            quant_dir / "SRR1_abundance.tsv"
        )

    def test_returns_none_on_empty_dir(self, tmp_path: Path) -> None:
        """An empty quantification directory yields None."""
        quant_dir = tmp_path / "SRR1"
        quant_dir.mkdir()
        assert find_quantification_file(quant_dir, sample_id="SRR1") is None
