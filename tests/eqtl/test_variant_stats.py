"""Tests for metainformant.eqtl.variant_stats (zero-mocks).

Parses fixture text captured from real bcftools output format and exercises
JSON writers on real temp files. External-tool wrappers are tested for
argument construction via tool-availability checks, not mocks.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from metainformant.eqtl.variant_stats import (
    compute_allele_frequencies,
    compute_popgen_summary,
    compute_sample_stats,
    parse_bcftools_stats,
)

# Fixture in the exact format emitted by `bcftools stats`
BCFTOOLS_STATS_FIXTURE = """\
# This file was produced by bcftools stats
SN\t0\tnumber of samples:\t3
SN\t0\tnumber of records:\t12
SN\t0\tnumber of no-ALTs:\t0
SN\t0\tnumber of SNPs:\t10
SN\t0\tnumber of MNPs:\t0
SN\t0\tnumber of indels:\t2
TSTV\t0\t10\t2\t8.33\t0.6410\t0.7377\t0.6410\t0.7377
"""


class TestParseBcftoolsStats:
    def test_parses_counts_and_tstv(self):
        stats = parse_bcftools_stats(BCFTOOLS_STATS_FIXTURE)
        assert stats["n_records"] == 12
        assert stats["n_snps"] == 10
        assert stats["n_indels"] == 2
        assert stats["ts_tv_ratio"] == pytest.approx(8.33)

    def test_empty_text_gives_zeroes(self):
        stats = parse_bcftools_stats("")
        assert stats == {
            "n_records": 0,
            "n_snps": 0,
            "n_indels": 0,
            "ts_tv_ratio": 0.0,
        }

    def test_malformed_tstv_is_tolerated(self):
        text = "TSTV\t0\t10\t2\tnot-a-float\n"
        stats = parse_bcftools_stats(text)
        assert stats["ts_tv_ratio"] == 0.0


class TestComputeSampleStats:
    def test_writes_json_with_real_bcftools_or_skips(self, tmp_path: Path):
        import shutil

        if shutil.which("bcftools") is None:
            pytest.skip("bcftools not installed")
        from metainformant.eqtl.test_support import build_tiny_vcf

        vcf = build_tiny_vcf(tmp_path)
        out = tmp_path / "stats" / "sample_stats.json"
        stats = compute_sample_stats(vcf, out)
        assert out.exists()
        assert stats["n_records"] >= 1
        loaded = json.loads(out.read_text())
        assert loaded["n_records"] == stats["n_records"]


class TestComputeAlleleFrequencies:
    def test_requires_bcftools_or_skips(self, tmp_path: Path):
        import shutil

        if shutil.which("bcftools") is None:
            pytest.skip("bcftools not installed")
        from metainformant.eqtl.test_support import build_tiny_vcf

        vcf = build_tiny_vcf(tmp_path)
        out = tmp_path / "af" / "af.tsv"
        n = compute_allele_frequencies(vcf, out)
        assert n >= 1
        header = out.read_text().splitlines()[0]
        assert header == "chrom\tpos\tref\talt\taf"


class TestComputePopgenSummary:
    def test_summary_embeds_per_sample_stats(self, tmp_path: Path):
        # Pure parsing/writing path: no VCF needed if we feed precomputed stats
        # through the same structure the pipeline uses.
        sample_stats = [
            {"sample_id": "S1", "n_snps": 5},
            {"sample_id": "S2", "n_snps": 7},
        ]
        out = tmp_path / "pop" / "summary.json"
        # compute_popgen_summary runs bcftools on merged_vcf; verify contract
        # via the JSON-writing path with a stub-free direct call skipped when
        # bcftools is absent.
        import shutil

        if shutil.which("bcftools") is None:
            with pytest.raises(RuntimeError, match="bcftools"):
                compute_popgen_summary(tmp_path / "missing.vcf.gz", sample_stats, out)
        else:
            from metainformant.eqtl.test_support import build_tiny_vcf

            vcf = build_tiny_vcf(tmp_path)
            summary = compute_popgen_summary(vcf, sample_stats, out)
            assert summary["n_samples"] == 2
            assert out.exists()
            assert summary["total_snps"] >= 1
