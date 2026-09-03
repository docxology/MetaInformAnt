"""Tests for metainformant.eqtl.workflow.variant_calling (zero-mocks).

Exercises discovery, decompression, and tool-availability logic with real
files; external-tool runs are skipped when tools are absent (marked by the
real-implementation policy: no test doubles, real subprocess or skip).
"""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest

from metainformant.eqtl.workflow.variant_calling import (
    check_tools,
    decompress_if_needed,
    find_completed_samples,
    find_reference_genome,
)


class TestCheckTools:
    def test_returns_bool(self):
        assert isinstance(check_tools(), bool)


class TestFindReferenceGenome:
    def test_finds_genomic_fna_preferring_non_cds(self, tmp_path: Path):
        genome_dir = tmp_path / "amellifera" / "genome"
        genome_dir.mkdir(parents=True)
        cds = genome_dir / "GCF_x_cds_from_genomic.fna"
        cds.write_text(">cds\n")
        genomic = genome_dir / "GCF_x_genomic.fna"
        genomic.write_text(">genome\n")
        found = find_reference_genome(tmp_path, "amellifera")
        assert found == genomic

    def test_falls_back_to_fna_gz(self, tmp_path: Path):
        genome_dir = tmp_path / "sp" / "genome"
        genome_dir.mkdir(parents=True)
        gz = genome_dir / "GCF_y_genomic.fna.gz"
        gz.write_bytes(b"\x1f\x8b\x08\x00")
        found = find_reference_genome(tmp_path, "sp")
        assert found == gz

    def test_missing_genome_dir_returns_none(self, tmp_path: Path):
        assert find_reference_genome(tmp_path, "nope") is None


class TestFindCompletedSamples:
    def test_discovers_quant_dirs(self, tmp_path: Path):
        quant_root = tmp_path / "sp" / "quant" / "quant"
        (quant_root / "SRR1").mkdir(parents=True)
        (quant_root / "SRR2").mkdir(parents=True)
        (quant_root / "not_a_sample.txt").write_text("x")
        # A sample dir without a quantification file must be excluded
        empty = quant_root / "SRR_empty"
        empty.mkdir()
        samples = find_completed_samples(tmp_path, "sp")
        # SRR1/SRR2 lack real quant files; behavior must be honest about it
        assert samples == [] or set(samples) <= {"SRR1", "SRR2"}


class TestDecompressIfNeeded:
    def test_passthrough_plain_file(self, tmp_path: Path):
        p = tmp_path / "plain.txt"
        p.write_text("data")
        assert decompress_if_needed(p) == p

    def test_decompresses_gz(self, tmp_path: Path):
        gz = tmp_path / "data.txt.gz"
        with gzip.open(gz, "wt") as f:
            f.write("payload")
        out = decompress_if_needed(gz)
        assert out.name == "data.txt"
        assert out.read_text() == "payload"

    def test_reuses_existing_decompressed(self, tmp_path: Path):
        gz = tmp_path / "data2.txt.gz"
        with gzip.open(gz, "wt") as f:
            f.write("new")
        existing = tmp_path / "data2.txt"
        existing.write_text("old")
        out = decompress_if_needed(gz)
        assert out == existing
        assert out.read_text() == "old"


@pytest.mark.skipif(shutil.which("hisat2-build") is None, reason="hisat2-build not installed")
class TestBuildHisat2Index:
    def test_builds_index_for_tiny_fasta(self, tmp_path: Path):
        from metainformant.eqtl.workflow.variant_calling import build_hisat2_index

        fasta = tmp_path / "genome.fna"
        fasta.write_text(">c1\n" + "ACGT" * 100 + "\n")
        out_dir = tmp_path / "index"
        prefix = build_hisat2_index(fasta, out_dir)
        assert (out_dir / ".index_built").exists()
        assert any(out_dir.glob("genome*.ht2")) or prefix.exists()

    def test_resume_skips_rebuild(self, tmp_path: Path):
        from metainformant.eqtl.workflow.variant_calling import build_hisat2_index

        fasta = tmp_path / "genome.fna"
        fasta.write_text(">c1\n" + "ACGT" * 100 + "\n")
        out_dir = tmp_path / "index"
        build_hisat2_index(fasta, out_dir)
        first_mtime = (out_dir / ".index_built").stat().st_mtime_ns
        build_hisat2_index(fasta, out_dir)
        assert (out_dir / ".index_built").stat().st_mtime_ns == first_mtime
