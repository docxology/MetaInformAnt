"""Zero-mocks tests for gwas.reporting.audit and gwas.validation.output_validator.

Real files on disk, real hashes, real in-memory config dicts.
"""

from __future__ import annotations

import hashlib
import json

from metainformant.gwas.reporting import audit
from metainformant.gwas.validation import output_validator


class TestSha256File:
    def test_known_digest(self, tmp_path) -> None:
        f = tmp_path / "data.bin"
        f.write_bytes(b"metainformant")
        expected = hashlib.sha256(b"metainformant").hexdigest()
        assert audit.sha256_file(f) == expected

    def test_large_file_chunked(self, tmp_path) -> None:
        f = tmp_path / "big.bin"
        payload = b"x" * (3 * 1024 * 1024)  # > default 1MB chunk
        f.write_bytes(payload)
        assert audit.sha256_file(f) == hashlib.sha256(payload).hexdigest()

    def test_small_chunk_size(self, tmp_path) -> None:
        f = tmp_path / "small.bin"
        f.write_bytes(b"abcdef")
        assert audit.sha256_file(f, chunk_size=2) == hashlib.sha256(b"abcdef").hexdigest()


class TestChecksumDirectory:
    def test_only_matching_extensions(self, tmp_path) -> None:
        (tmp_path / "plot.png").write_bytes(b"png")
        (tmp_path / "results.tsv").write_text("col\n")
        (tmp_path / "notes.txt").write_text("skip")
        (tmp_path / "sub").mkdir()
        (tmp_path / "sub" / "nested.json").write_text("{}")
        sums = audit.checksum_directory(tmp_path)
        assert set(sums.keys()) == {"plot.png", "results.tsv", str("sub/nested.json")}
        assert sums["plot.png"] == hashlib.sha256(b"png").hexdigest()

    def test_custom_extensions(self, tmp_path) -> None:
        (tmp_path / "a.vcf").write_text("##fileformat\n")
        (tmp_path / "b.png").write_bytes(b"png")
        sums = audit.checksum_directory(tmp_path, extensions=[".vcf"])
        assert set(sums.keys()) == {"a.vcf"}


class TestLibraryVersions:
    def test_reports_installed(self) -> None:
        versions = audit.library_versions(["numpy", "metainformant"])
        assert versions["numpy"] != "not installed"
        assert versions["metainformant"] != "not installed"

    def test_reports_missing(self) -> None:
        versions = audit.library_versions(["definitely_not_a_real_lib_xyz"])
        assert versions["definitely_not_a_real_lib_xyz"] == "not installed"


class TestGenerateAuditLog:
    def test_audit_log_structure(self, tmp_path) -> None:
        config = {
            "gwas": {"model": "mixed", "significance_threshold": 5e-8, "n_pca_components": 10},
            "qc": {"min_maf": 0.02, "max_missing_rate": 0.1, "hwe_p_threshold": 1e-6},
            "ld_pruning": {"r2_threshold": 0.2},
            "post_gwas": {"top_hits_to_label": 5, "gene_annotation_radius_bp": 50000},
            "samples": {"phenotypes": {"default_trait": "wing"}},
        }
        (tmp_path / "res").mkdir()
        (tmp_path / "res" / "out.tsv").write_text("x\t1\n")
        (tmp_path / "proc").mkdir()
        result = audit.generate_audit_log(
            config=config,
            results_dir=tmp_path / "res",
            processed_dir=tmp_path / "proc",
            raw_dirs=[tmp_path],
            variant_count=1234,
            sig_count=5,
            lambda_gc=1.04,
        )
        assert isinstance(result, dict)
        assert result["pipeline_version"] == "metainformant-gwas"
        assert result["parameters"]["trait"] == "wing"
        assert result["parameters"]["model"] == "mixed"
        assert result["parameters"]["significance_threshold"] == 5e-8
        assert "environment" in result and "run_timestamp_utc" in result


class TestValidatePipelineOutputs:
    @staticmethod
    def _write_outputs(results_dir, *, lambda_gc: float = 1.05, n_gw: int = 3, sig_p: float = 1e-9) -> None:
        trait_dir = results_dir / "wing" / "mixed"
        post = trait_dir / "post_gwas"
        post.mkdir(parents=True)
        headers = ["CHR", "POS", "SNP", "BETA", "SE", "P", "MAF"]
        summary_rows = [("1", 100, "s1", 0.1, 0.01, sig_p, 0.3), ("1", 200, "s2", 0.02, 0.01, 0.5, 0.2)]
        with open(trait_dir / "summary_statistics.tsv", "w") as f:
            f.write("\t".join(headers) + "\n")
            for r in summary_rows:
                f.write("\t".join(map(str, r)) + "\n")
        with open(trait_dir / "significant_hits.tsv", "w") as f:
            f.write("\t".join(headers) + "\n")
            f.write("\t".join(map(str, summary_rows[0])) + "\n")
        for name in [
            "manhattan_plot.png",
            "qq_plot.png",
            "pca_plot.png",
            "kinship_plot.png",
            "effect_size_plot.png",
            "maf_spectrum_plot.png",
        ]:
            (trait_dir / name).write_bytes(b"png")
        (post / "volcano_plot.png").write_bytes(b"png")
        (post / "qq_stratified.png").write_bytes(b"png")
        (post / "chrom_summary.tsv").write_text("chrom\tn\n")
        (post / "post_gwas_results.json").write_text(
            json.dumps(
                {
                    "summary_statistics": {
                        "p_value_calibration": {"lambda_gc": lambda_gc},
                        "significance_counts": {"genome_wide_5e-8": n_gw},
                    },
                    "n_variants_analyzed": 2,
                    "sign_test": {"p": 0.4},
                    "credible_sets": [],
                }
            )
        )

    @staticmethod
    def _config(results_dir) -> dict:
        return {
            "paths": {"results_dir": str(results_dir), "processed_dir": "x", "phenotype_dir": "y", "raw_dir": "z"},
            "samples": {"phenotypes": {"default_trait": "wing"}},
            "gwas": {"model": "mixed", "significance_threshold": 5e-8},
        }

    def test_all_valid_outputs_pass(self, tmp_path) -> None:
        self._write_outputs(tmp_path)
        report = output_validator.validate_pipeline_outputs(self._config(tmp_path))
        assert report["all_ok"] is True
        assert report["n_fail"] == 0
        labels = {c["label"] for c in report["checks"]}
        assert "summary_statistics.tsv schema" in labels
        assert "significant_hits.tsv consistency" in labels
        assert "lambda_gc finite and positive" in labels
        assert "post_gwas_results.json schema" in labels

    def test_missing_file_fails(self, tmp_path) -> None:
        self._write_outputs(tmp_path)
        (tmp_path / "wing" / "mixed" / "manhattan_plot.png").unlink()
        report = output_validator.validate_pipeline_outputs(self._config(tmp_path))
        assert report["all_ok"] is False
        assert report["n_fail"] >= 1

    def test_inflated_lambda_warns_not_fails(self, tmp_path) -> None:
        self._write_outputs(tmp_path, lambda_gc=3.5)
        report = output_validator.validate_pipeline_outputs(self._config(tmp_path))
        lam = [c for c in report["checks"] if c["label"] == "lambda_gc finite and positive"][0]
        assert lam["status"] == "WARN"
        assert report["all_ok"] is True  # warnings do not fail the gate

    def test_hit_above_threshold_fails_consistency(self, tmp_path) -> None:
        self._write_outputs(tmp_path, sig_p=0.01)  # not genome-wide significant
        report = output_validator.validate_pipeline_outputs(self._config(tmp_path))
        consistency = [c for c in report["checks"] if c["label"] == "significant_hits.tsv consistency"][0]
        assert consistency["status"] == "FAIL"

    def test_require_hits_flag(self, tmp_path) -> None:
        self._write_outputs(tmp_path, n_gw=0)
        relaxed = output_validator.validate_pipeline_outputs(self._config(tmp_path), require_hits=False)
        assert relaxed["all_ok"] is True
        strict = output_validator.validate_pipeline_outputs(self._config(tmp_path), require_hits=True)
        assert strict["all_ok"] is False

    def test_missing_outputs_entirely(self, tmp_path) -> None:
        report = output_validator.validate_pipeline_outputs(self._config(tmp_path))
        assert report["all_ok"] is False
        assert report["n_pass"] == 0
