"""Zero-mocks tests for gwas.data.vcf_utils (discovery, counting, IDs)."""

from __future__ import annotations

import shutil
import subprocess

import pytest

from metainformant.gwas.data import vcf_utils

BCFTOOLS = shutil.which("bcftools")


MINIMAL_VCF = """##fileformat=VCFv4.2
##contig=<ID=NC_037638.1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
NC_037638.1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1
NC_037638.1\t200\t.\tC\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/1
NC_037638.1\t300\t.\tT\tA\t50\tPASS\t.\tGT\t1/1\t0/1\t0/0
NC_037639.1\t100\t.\tG\tC\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0
"""


class TestDiscoverSampleVcfs:
    def test_discovers_and_excludes_artifacts(self, tmp_path) -> None:
        for name in ("S1.vcf", "S2.vcf", "S1_merged.vcf", "subset_S3.vcf"):
            (tmp_path / name).write_text("##fileformat=VCFv4.2\n")
        found = vcf_utils.discover_sample_vcfs(tmp_path)
        assert [p.name for p in found] == ["S1.vcf", "S2.vcf"]

    def test_empty_dir(self, tmp_path) -> None:
        assert vcf_utils.discover_sample_vcfs(tmp_path) == []


class TestCountVariants:
    def test_counts_data_lines_only(self, tmp_path) -> None:
        if not shutil.which("bcftools"):
            pytest.skip("bcftools not installed")
        f = tmp_path / "in.vcf"
        f.write_text(MINIMAL_VCF)
        assert vcf_utils.count_variants(f) == 4

    def test_counts_gzipped(self, tmp_path) -> None:
        if not BCFTOOLS:
            pytest.skip("bcftools not installed")
        plain = tmp_path / "in.vcf"
        plain.write_text(MINIMAL_VCF)
        out = vcf_utils.bgzip_and_index(plain)
        assert out is not None and out.exists()
        assert vcf_utils.count_variants(out) == 4


class TestBgzipAndIndex:
    def test_creates_bgz_and_tbi(self, tmp_path) -> None:
        if not BCFTOOLS:
            pytest.skip("bcftools not installed")
        plain = tmp_path / "in.vcf"
        plain.write_text(MINIMAL_VCF)
        out = vcf_utils.bgzip_and_index(plain)
        assert out is not None
        assert out.name == "in.vcf.gz"
        assert (tmp_path / "in.vcf.gz.tbi").exists()
        # Content preserved
        check = subprocess.run([BCFTOOLS, "view", "-H", str(out)], capture_output=True, text=True)
        assert len(check.stdout.strip().splitlines()) == 4


class TestMergeVcfs:
    def test_merges_two_vcfs(self, tmp_path) -> None:
        if not BCFTOOLS:
            pytest.skip("bcftools not installed")
        inputs = []
        for name in ("a", "b"):
            plain = tmp_path / f"{name}.vcf"
            plain.write_text(MINIMAL_VCF)
            gz = vcf_utils.bgzip_and_index(plain)
            assert gz is not None
            inputs.append(gz)
        out = tmp_path / "merged.vcf.gz"
        result = vcf_utils.merge_vcfs(inputs, out)
        assert result is not None and out.exists()
        check = subprocess.run([BCFTOOLS, "view", "-H", str(out)], capture_output=True, text=True)
        assert check.returncode == 0
        assert len(check.stdout.strip().splitlines()) > 0


class TestSubsampleVcf:
    def test_subsampled_variant_count(self, tmp_path) -> None:
        if not BCFTOOLS:
            pytest.skip("bcftools not installed")
        plain = tmp_path / "in.vcf"
        plain.write_text(MINIMAL_VCF)
        out = tmp_path / "sub.vcf"
        result = vcf_utils.subsample_vcf(plain, out, fraction=0.5, seed=42)
        assert result is not None and out.exists()
        n = vcf_utils.count_variants(out)
        assert 0 <= n <= 4


class TestExtractSampleIds:
    def test_extracts_header_samples(self, tmp_path) -> None:
        if not BCFTOOLS:
            pytest.skip("bcftools not installed")
        f = tmp_path / "in.vcf"
        f.write_text(MINIMAL_VCF)
        ids = vcf_utils.extract_sample_ids(f)
        assert ids == ["S1", "S2", "S3"]
