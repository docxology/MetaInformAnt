"""Regression tests: eQTL code that reported a placeholder or silently dropped input now reports what happened.

* ``parse_bcftools_stats`` returned ``ts_tv_ratio = 0.0`` when bcftools printed no ratio, or ``nan`` (it prints ``nan``
  when there is no transversion), and ``run_pipeline`` copied that into its summary; a ratio of 0 means "no
  transitions", which is not what any of these cases says. It is now None, and so are the SNP totals of a run in
  which no sample completed;
* ``create_synthetic_data`` validated ``n_variants >= 5 * n_genes`` but always generated ``5 * n_genes`` variants, so
  asking for more silently gave fewer rows than asked;
* ``align_reads`` aligned a three-file run (``<run>_1``, ``<run>_2`` and the unpaired ``<run>``) as single-end using the
  first file only, which dropped the read pair without a message;
* ``merge_vcfs`` failed on a missing destination directory when given a single VCF (a one-sample run), because the
  directory was only created for two or more inputs.

Oracles are hand-built text in the exact format bcftools emits, and hand-counted values.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from metainformant.eqtl.variant_stats import parse_bcftools_stats
from metainformant.eqtl.workflow import variant_calling
from metainformant.eqtl.workflow.synthetic import create_synthetic_data

# ---------------------------------------------------------------------------
# Ti/Tv: absent is not zero
# ---------------------------------------------------------------------------


def tstv(ratio: str) -> str:
    return f"SN\t0\tnumber of SNPs:\t4\nTSTV\t0\t4\t0\t{ratio}\t0.5\t0.0\t{ratio}\n"


@pytest.mark.parametrize("printed", ["nan", "-nan", "NaN", "inf", "-inf", "not-a-number", ""])
def test_a_ratio_that_is_not_a_finite_number_is_none_not_zero(printed: str) -> None:
    assert parse_bcftools_stats(tstv(printed))["ts_tv_ratio"] is None


@pytest.mark.parametrize("printed, expected", [("2.15", 2.15), ("0.00", 0.0), ("10", 10.0)])
def test_a_ratio_bcftools_reports_is_kept_including_a_real_zero(printed: str, expected: float) -> None:
    ratio = parse_bcftools_stats(tstv(printed))["ts_tv_ratio"]

    assert ratio == pytest.approx(expected)
    assert ratio is not None


def test_no_tstv_line_is_no_ratio() -> None:
    stats = parse_bcftools_stats("SN\t0\tnumber of records:\t3\nSN\t0\tnumber of SNPs:\t3\n")

    assert stats["n_records"] == 3 and stats["n_snps"] == 3
    assert stats["ts_tv_ratio"] is None


def test_a_later_valid_tstv_line_wins_over_an_earlier_invalid_one() -> None:
    text = "TSTV\t0\t4\t0\tnan\t0\t0\tnan\nTSTV\t1\t8\t4\t2.00\t0\t0\t2.00\n"

    assert parse_bcftools_stats(text)["ts_tv_ratio"] == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# create_synthetic_data: n_variants is the number of variants
# ---------------------------------------------------------------------------


class TestSyntheticVariantCount:
    def test_the_requested_number_of_variants_is_generated(self) -> None:
        _, genotypes, _, variants = create_synthetic_data(n_genes=4, n_variants=30, n_samples=12, seed=7)

        assert genotypes.shape == (30, 12)
        assert len(variants) == 30
        assert variants["variant_id"].is_unique
        assert list(genotypes.index) == list(variants["variant_id"])

    def test_the_cis_variants_are_the_same_as_without_background_variants(self) -> None:
        _, plain, _, plain_variants = create_synthetic_data(n_genes=4, n_variants=20, n_samples=12, seed=7)
        _, extended, _, extended_variants = create_synthetic_data(n_genes=4, n_variants=30, n_samples=12, seed=7)

        pd.testing.assert_frame_equal(extended.iloc[:20], plain)
        pd.testing.assert_frame_equal(extended_variants.iloc[:20].reset_index(drop=True), plain_variants)

    def test_background_variants_are_outside_every_cis_window(self) -> None:
        _, _, genes, variants = create_synthetic_data(n_genes=4, n_variants=30, n_samples=12, seed=7)
        background = variants.iloc[20:]

        assert len(background) == 10
        window = 1_000_000
        for _, gene in genes.iterrows():
            same_chromosome = background[background["chrom"] == gene["chrom"]]
            assert (abs(same_chromosome["position"] - gene["tss_position"]) > window).all()
        assert set(background["chrom"]) == {"1", "2", "3", "4", "5"}

    def test_background_variants_do_not_drive_any_gene(self) -> None:
        expression, genotypes, _, _ = create_synthetic_data(n_genes=30, n_variants=200, n_samples=80, seed=5)

        strongest = max(
            abs(np.corrcoef(expression.iloc[g].to_numpy(), genotypes.iloc[row].to_numpy())[0, 1])
            for g in range(30)
            for row in range(150, 200)
        )
        assert strongest < 0.5  # a true cis effect in this setup correlates at 0.3-0.9; chance does not reach 0.5

    def test_the_exact_count_still_works_and_too_few_still_raises(self) -> None:
        _, genotypes, _, variants = create_synthetic_data(n_genes=3, n_variants=15, n_samples=10)
        assert genotypes.shape[0] == 15 and len(variants) == 15

        with pytest.raises(ValueError, match="n_variants"):
            create_synthetic_data(n_genes=3, n_variants=14, n_samples=10)


# ---------------------------------------------------------------------------
# hisat2 command line
# ---------------------------------------------------------------------------

INDEX = Path("/idx/genome")
BASE = ["hisat2", "--dta", "-p", "4", "-x", "/idx/genome"]


class TestHisat2Command:
    def test_one_file_is_single_end(self) -> None:
        assert variant_calling.hisat2_command([Path("SRR1.fastq.gz")], INDEX) == BASE + ["-U", "SRR1.fastq.gz"]

    def test_two_files_are_a_read_pair(self) -> None:
        files = [Path("SRR1_1.fastq.gz"), Path("SRR1_2.fastq.gz")]

        assert variant_calling.hisat2_command(files, INDEX) == BASE + ["-1", "SRR1_1.fastq.gz", "-2", "SRR1_2.fastq.gz"]

    def test_a_pair_plus_the_unpaired_reads_uses_all_three_files(self) -> None:
        """Sorted as the downloader returns them, the unpaired file comes first: it used to be the only one used."""
        files = sorted([Path("SRR1_2.fastq.gz"), Path("SRR1.fastq.gz"), Path("SRR1_1.fastq.gz")])
        assert files[0].name == "SRR1.fastq.gz"

        command = variant_calling.hisat2_command(files, INDEX, threads=8)

        assert command == [
            "hisat2",
            "--dta",
            "-p",
            "8",
            "-x",
            "/idx/genome",
            "-1",
            "SRR1_1.fastq.gz",
            "-2",
            "SRR1_2.fastq.gz",
            "-U",
            "SRR1.fastq.gz",
        ]

    def test_several_unpaired_files_are_joined_for_the_unpaired_option(self) -> None:
        files = [Path("a.fastq.gz"), Path("SRR1_1.fq.gz"), Path("SRR1_2.fq.gz"), Path("b.fastq.gz")]

        command = variant_calling.hisat2_command(files, INDEX)

        assert command[-2:] == ["-U", "a.fastq.gz,b.fastq.gz"]
        assert command[6:10] == ["-1", "SRR1_1.fq.gz", "-2", "SRR1_2.fq.gz"]

    @pytest.mark.parametrize(
        "names",
        [
            ["a.fastq.gz", "b.fastq.gz", "c.fastq.gz"],  # no mates at all
            ["SRR1_1.fastq.gz", "SRR1_1.fastq.gz", "SRR1_2.fastq.gz"],  # two first mates
            ["SRR1_1.fastq.gz", "SRR2_1.fastq.gz", "SRR1.fastq.gz"],  # no second mate
        ],
    )
    def test_an_ambiguous_layout_is_an_error_not_a_guess_that_drops_reads(self, names: list) -> None:
        with pytest.raises(ValueError, match="cannot tell the read layout of 3 FASTQ files"):
            variant_calling.hisat2_command([Path(n) for n in names], INDEX)

    def test_no_files_is_an_error(self) -> None:
        with pytest.raises(ValueError, match="no FASTQ files"):
            variant_calling.hisat2_command([], INDEX)


# ---------------------------------------------------------------------------
# merge_vcfs with one sample
# ---------------------------------------------------------------------------


def test_a_single_vcf_is_copied_into_a_destination_directory_that_does_not_exist_yet(tmp_path: Path) -> None:
    """The copy comes first; indexing then fails (no bcftools, or not a real VCF), which is all this test needs."""
    source = tmp_path / "sample.vcf.gz"
    source.write_bytes(b"not a real vcf")
    destination = tmp_path / "population" / "merged.vcf.gz"

    with pytest.raises((FileNotFoundError, subprocess.CalledProcessError)):
        variant_calling.merge_vcfs([source], destination)

    assert destination.read_bytes() == b"not a real vcf"  # it used to fail on the missing directory before copying
    assert shutil.which("bcftools") is None or destination.parent.is_dir()


def test_no_vcf_is_not_a_merge(tmp_path: Path) -> None:
    assert variant_calling.merge_vcfs([], tmp_path / "population" / "merged.vcf.gz") is False
