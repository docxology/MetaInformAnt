"""Zero-mocks tests for gwas.simulation.synthetic (analytical power + synthetic GWAS data)."""

from __future__ import annotations

import json
import math
import subprocess

import pytest

from metainformant.gwas.simulation import synthetic

BCFTOOLS = None
try:
    from shutil import which

    BCFTOOLS = which("bcftools")
except ImportError:
    pass


class TestNormalCdfInverse:
    def test_inverse_pairing(self) -> None:
        for p in [0.01, 0.1, 0.5, 0.9, 0.975]:
            z = synthetic.inv_normal_cdf(p)
            assert abs(synthetic.normal_cdf(z) - p) < 1e-4

    def test_known_quantiles(self) -> None:
        assert abs(synthetic.inv_normal_cdf(0.975) - 1.959964) < 1e-3
        assert abs(synthetic.inv_normal_cdf(0.5)) < 1e-6

    def test_out_of_range_raises(self) -> None:
        with pytest.raises(ValueError):
            synthetic.inv_normal_cdf(0.0)
        with pytest.raises(ValueError):
            synthetic.inv_normal_cdf(1.0)

    def test_symmetry(self) -> None:
        assert abs(synthetic.inv_normal_cdf(0.3) + synthetic.inv_normal_cdf(0.7)) < 1e-9


class TestCalculatePower:
    def test_matches_manual_ncp_formula(self) -> None:
        beta, maf, n = 0.4, 0.25, 1000
        ncp = abs(beta) * math.sqrt(n * 2 * maf * (1 - maf))
        expected = 1.0 - synthetic.normal_cdf(synthetic.inv_normal_cdf(1 - 5e-8 / 2) - ncp)
        assert abs(synthetic.calculate_power(beta, maf, n) - expected) < 1e-9

    def test_zero_effect_zero_power(self) -> None:
        assert synthetic.calculate_power(0.0, 0.5, 5000) < 1e-6

    def test_power_increases_with_beta_and_n(self) -> None:
        p1 = synthetic.calculate_power(0.1, 0.3, 500)
        p2 = synthetic.calculate_power(0.5, 0.3, 500)
        p3 = synthetic.calculate_power(0.5, 0.3, 5000)
        assert p2 > p1
        assert p3 > p2

    def test_power_bounded(self) -> None:
        p = synthetic.calculate_power(3.0, 0.5, 100000)
        assert 0.0 <= p <= 1.0


class TestCalibrationModes:
    def test_modes_documented(self) -> None:
        assert set(synthetic.CALIBRATION_MODES) == {"NULL", "LOW", "HIGH"}
        assert synthetic.CALIBRATION_MODES["NULL"]["h2"] == 0.0
        assert synthetic.CALIBRATION_MODES["HIGH"]["h2"] > synthetic.CALIBRATION_MODES["LOW"]["h2"]


class TestCreateSyntheticData:
    def test_generates_vcf_phenotype_and_params(self, tmp_path) -> None:
        params = synthetic.create_synthetic_data(
            tmp_path, n_samples=12, n_variants=100, calibration_mode="LOW", seed=42
        )
        vcf = tmp_path / "synthetic.vcf"
        pheno = tmp_path / "phenotype.csv"
        meta = tmp_path / "synthetic_params.json"
        assert vcf.exists() and pheno.exists() and meta.exists()
        assert params["n_samples"] == 12
        assert params["n_variants"] == 100
        assert params["n_causal"] >= 1
        with open(meta) as f:
            loaded = json.load(f)
        assert loaded["seed"] == 42
        assert loaded["calibration_mode"] == "LOW"

    def test_vcf_parseable_and_correct_dimensions(self, tmp_path) -> None:
        synthetic.create_synthetic_data(tmp_path, n_samples=8, n_variants=50, calibration_mode="NULL", seed=7)
        vcf = tmp_path / "synthetic.vcf"
        header_samples = []
        n_variants = 0
        with open(vcf) as f:
            for line in f:
                if line.startswith("#CHROM"):
                    header_samples = line.strip().split("\t")[9:]
                elif not line.startswith("#"):
                    n_variants += 1
                    fields = line.strip().split("\t")
                    assert len(fields) == 9 + 8  # fixed cols + 8 samples
        assert header_samples == [f"C{i}WORK" for i in range(8)] or len(header_samples) == 8
        assert n_variants == 50

    def test_null_mode_has_zero_effects(self, tmp_path) -> None:
        params = synthetic.create_synthetic_data(tmp_path, n_samples=10, n_variants=60, calibration_mode="NULL", seed=3)
        assert params["h2"] == 0.0
        # NULL: all phenotypes should be pure noise, params record zero causal power basis
        assert params["n_causal"] == 0 or params["mean_causal_power"] == 0.0

    def test_deterministic_with_same_seed(self, tmp_path) -> None:
        d1 = tmp_path / "a"
        d2 = tmp_path / "b"
        synthetic.create_synthetic_data(d1, n_samples=10, n_variants=80, seed=99)
        synthetic.create_synthetic_data(d2, n_samples=10, n_variants=80, seed=99)
        assert (d1 / "synthetic.vcf").read_bytes() == (d2 / "synthetic.vcf").read_bytes()
        assert (d1 / "phenotype.csv").read_bytes() == (d2 / "phenotype.csv").read_bytes()

    def test_explicit_sample_ids_preserved(self, tmp_path) -> None:
        ids = ["C15G", "C16G", "I12ITQ", "M03WORK", "R21G", "R22G"]
        synthetic.create_synthetic_data(tmp_path, sample_ids=ids, n_variants=30, seed=5)
        content = (tmp_path / "synthetic.vcf").read_text()
        header = [ln for ln in content.splitlines() if ln.startswith("#CHROM")][0]
        assert header.split("\t")[9:] == ids

    def test_bcftools_can_read_output(self, tmp_path) -> None:
        if not BCFTOOLS:
            pytest.skip("bcftools not installed")
        synthetic.create_synthetic_data(tmp_path, n_samples=6, n_variants=40, seed=11)
        vcf_path = str(tmp_path / "synthetic.vcf")
        check = subprocess.run([BCFTOOLS, "view", "-H", vcf_path], capture_output=True, text=True)
        assert check.returncode == 0
        assert len(check.stdout.strip().splitlines()) == 40
