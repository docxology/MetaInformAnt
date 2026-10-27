"""Real-data tests for the P. barbatus pipeline entry points.

``pbarbatus_comprehensive``, ``pbarbatus_full_scale`` and ``pbarbatus_end_to_end`` each expose
``run(config_path)``. They are the business logic behind ``scripts/gwas/pipelines/*.py``. These tests
drive every ``run`` on a small, real, synthetic cohort (VCF + phenotype table + YAML config written
to ``tmp_path``) and check the files the pipelines write.

Contracts pinned here (all three modules unless noted):

* ``run(config_path)`` uses the config it is given, not the repository default.
* QC goes through ``apply_qc_filters(vcf_data, qc_config)`` and every later step uses the filtered
  variants and genotypes. Thresholds that YAML delivers as strings (``hwe_pval: 1e-6`` loads as the
  string ``"1e-6"``) are coerced to floats, and ``work_dir`` strings become ``Path`` objects before
  they are joined.
* Association rows carry the VCF record fields under their parsed lower-case names
  (``id``/``chrom``/``pos``), so the tables and Manhattan plots have coordinates.
* ``pbarbatus_full_scale``: ``pca_plot``/``kinship_heatmap``/``maf_distribution`` are called with
  their real signatures, the MAF comes from the non-missing dosages, and the closing file listing
  copes with the three-field trait tuples.
* ``pbarbatus_end_to_end``: significance flags come from the corrections' per-test booleans,
  Bonferroni-adjusted p-values stay within [0, 1], and the scree plot receives an array.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Dict, List

import numpy as np
import pandas as pd
import pytest

from metainformant.gwas.workflow import pbarbatus_comprehensive, pbarbatus_end_to_end, pbarbatus_full_scale

N_SAMPLES = 60
N_POLYMORPHIC = 80
N_MONOMORPHIC = 6
CAUSAL_INDEX = 10  # rs10 = chr1:11000 drives body_size_mm and disease_resistance
SEVEN_TRAITS = [
    "body_size_mm",
    "foraging_activity",
    "colony_size",
    "temp_tolerance_C",
    "longevity_days",
    "disease_resistance",
    "caste",
]


@dataclass(frozen=True)
class PipelineSpec:
    """Where a pipeline writes its outputs, relative to the config's ``work_dir``."""

    module: ModuleType
    results_subdir: str
    plots_subdir: str
    assoc_template: str
    traits: List[str]
    plot_files: List[str]


SPECS: Dict[str, PipelineSpec] = {
    "comprehensive": PipelineSpec(
        module=pbarbatus_comprehensive,
        results_subdir="results",
        plots_subdir="plots",
        assoc_template="{trait}_associations.tsv",
        traits=["body_size_mm", "foraging_activity", "disease_resistance"],
        plot_files=[
            f"{kind}_{trait}.png"
            for trait in ["body_size_mm", "foraging_activity", "disease_resistance"]
            for kind in ("manhattan", "qq")
        ],
    ),
    "full_scale": PipelineSpec(
        module=pbarbatus_full_scale,
        results_subdir="results_full",
        plots_subdir="plots_full",
        assoc_template="{trait}_associations_full.tsv",
        traits=SEVEN_TRAITS,
        plot_files=[f"{kind}_{trait}_full.png" for trait in SEVEN_TRAITS for kind in ("manhattan", "qq")]
        + ["pca_scatter_full.png", "kinship_heatmap_full.png", "maf_distribution_full.png"],
    ),
    "end_to_end": PipelineSpec(
        module=pbarbatus_end_to_end,
        results_subdir="results",
        plots_subdir="plots",
        assoc_template="{trait}_associations.tsv",
        traits=SEVEN_TRAITS,
        plot_files=[f"{kind}_{trait}.png" for trait in SEVEN_TRAITS for kind in ("manhattan", "qq")]
        + [
            "kinship_heatmap.png",
            "pca_scree_plot.png",
            "pca_scatter_PC1_vs_PC2.png",
            "pca_scatter_PC1_vs_PC3.png",
            "pca_scatter_PC2_vs_PC3.png",
        ],
    ),
}


def build_cohort(root: Path, n_polymorphic: int = N_POLYMORPHIC, n_monomorphic: int = N_MONOMORPHIC) -> Path:
    """Write a synthetic cohort and its YAML config below ``root``; return the config path.

    Polymorphic variants are sampled in Hardy-Weinberg proportions. ``rs{CAUSAL_INDEX}`` carries a
    strong additive effect on ``body_size_mm`` and a logistic effect on ``disease_resistance``. The
    ``mono*`` variants are fixed for the reference allele, so the MAF filter must remove them.
    """
    root.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(11)
    samples = [f"PB_{i:03d}" for i in range(N_SAMPLES)]
    freqs = rng.uniform(0.2, 0.5, max(n_polymorphic, 1))
    dosage = np.stack([rng.binomial(2, p, N_SAMPLES) for p in freqs[:n_polymorphic]]) if n_polymorphic else None

    causal = dosage[CAUSAL_INDEX] if dosage is not None else rng.binomial(2, 0.3, N_SAMPLES)
    body_size = 5.0 + 1.5 * causal + rng.normal(0.0, 0.4, N_SAMPLES)
    noise = rng.normal(0.0, 1.0, (4, N_SAMPLES))
    logit = -0.5 + 2.0 * (causal - causal.mean())
    disease = (rng.random(N_SAMPLES) < 1.0 / (1.0 + np.exp(-logit))).astype(int)
    caste = rng.integers(0, 2, N_SAMPLES)

    genotype_text = {0: "0/0", 1: "0/1", 2: "1/1"}
    vcf_path = root / "cohort.vcf"
    with open(vcf_path, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for v in range(n_polymorphic):
            chrom = "chr1" if v < n_polymorphic // 2 else "chr2"
            calls = "\t".join(genotype_text[int(g)] for g in dosage[v])
            handle.write(f"{chrom}\t{1000 * (v + 1)}\trs{v}\tA\tG\t60\tPASS\t.\tGT\t{calls}\n")
        for m in range(n_monomorphic):
            calls = "\t".join(["0/0"] * N_SAMPLES)
            handle.write(f"chr2\t{1_000_000 + 1000 * m}\tmono{m}\tA\tG\t60\tPASS\t.\tGT\t{calls}\n")

    pheno_path = root / "phenotypes.tsv"
    pheno = pd.DataFrame(
        {
            "sample_id": samples,
            "age_days": 30 + np.arange(N_SAMPLES),
            "colony": np.arange(N_SAMPLES) % 3,
            "sex": "F",
            "body_size_mm": body_size,
            "foraging_activity": noise[0],
            "colony_size": noise[1],
            "temp_tolerance_C": noise[2],
            "longevity_days": noise[3],
            "disease_resistance": disease,
            "caste": caste,
        }
    )
    pheno.to_csv(pheno_path, sep="\t", index=False)

    work_dir = root / "work"
    config_path = root / "config.yaml"
    config_path.write_text(f"""work_dir: {work_dir}
log_dir: {work_dir}/logs
threads: 1
genome:
  accession: GCF_TEST
  fasta: {root}/missing.fna
variants:
  vcf_files:
    - {vcf_path}
qc:
  min_maf: 0.05
  max_missing: 0.10
  hwe_pval: 1e-6
  min_qual: 30.0
samples:
  phenotype_file: {pheno_path}
structure:
  n_components: 5
  kinship_method: vanraden
output:
  results_dir: {work_dir}/results
  plots_dir: {work_dir}/plots
""")
    return config_path


@dataclass(frozen=True)
class PipelineRun:
    spec: PipelineSpec
    root: Path
    results_dir: Path
    plots_dir: Path

    def assoc(self, trait: str) -> pd.DataFrame:
        return pd.read_csv(self.results_dir / self.spec.assoc_template.format(trait=trait), sep="\t")


@pytest.fixture(scope="module")
def pipeline_runs(tmp_path_factory: pytest.TempPathFactory) -> Dict[str, PipelineRun]:
    """Run each pipeline once on its own copy of the cohort; the tests below inspect what it wrote."""
    runs: Dict[str, PipelineRun] = {}
    for name, spec in SPECS.items():
        root = tmp_path_factory.mktemp(f"pbarbatus_{name}")
        spec.module.run(str(build_cohort(root)))
        work_dir = root / "work"
        runs[name] = PipelineRun(spec, root, work_dir / spec.results_subdir, work_dir / spec.plots_subdir)
    return runs


@pytest.mark.parametrize("name", sorted(SPECS))
def test_pipeline_writes_outputs_next_to_the_supplied_config(name: str, pipeline_runs: Dict[str, PipelineRun]) -> None:
    """``run(config_path)`` must use that config: every table and plot lands under its work_dir."""
    run = pipeline_runs[name]
    for trait in run.spec.traits:
        table = run.results_dir / run.spec.assoc_template.format(trait=trait)
        assert table.stat().st_size > 0, table
    for plot_name in run.spec.plot_files:
        plot = run.plots_dir / plot_name
        assert plot.exists() and plot.stat().st_size > 0, plot
    for file_name in ("qc_results.json", "pca_results.json", "kinship_results.json"):
        assert (run.results_dir / file_name).stat().st_size > 0, file_name


@pytest.mark.parametrize("name", sorted(SPECS))
def test_pipeline_association_rows_carry_variant_coordinates(name: str, pipeline_runs: Dict[str, PipelineRun]) -> None:
    """Rows carry the VCF id/chrom/pos; the causal variant is the top body-size hit."""
    table = pipeline_runs[name].assoc("body_size_mm")

    assert table[["variant_id", "chrom", "pos"]].notna().all().all()
    top = table.sort_values("p_value").iloc[0]
    assert (top["variant_id"], top["chrom"], int(top["pos"])) == (
        f"rs{CAUSAL_INDEX}",
        "chr1",
        1000 * (CAUSAL_INDEX + 1),
    )
    assert top["beta"] == pytest.approx(1.5, abs=0.25)
    assert top["p_value"] < 1e-10


@pytest.mark.parametrize("name", sorted(SPECS))
def test_pipeline_tests_only_variants_that_pass_qc(name: str, pipeline_runs: Dict[str, PipelineRun]) -> None:
    """The monomorphic variants fail the MAF filter and must not reach association testing."""
    run = pipeline_runs[name]
    qc = json.loads((run.results_dir / "qc_results.json").read_text())
    assert qc["status"] == "success"
    assert qc["num_variants_before"] == N_POLYMORPHIC + N_MONOMORPHIC
    assert qc["num_variants_after"] == N_POLYMORPHIC
    assert qc["qc_filters"]["min_hwe_p"] == pytest.approx(1e-6)
    assert "filtered_data" not in qc  # the summary file is not a dump of the filtered genotypes

    for trait in run.spec.traits:
        table = run.assoc(trait)
        assert len(table) == N_POLYMORPHIC, trait
        assert not table["variant_id"].astype(str).str.startswith("mono").any(), trait


def test_end_to_end_significance_flags_and_adjusted_p_values(pipeline_runs: Dict[str, PipelineRun]) -> None:
    """Flags come from the corrections' per-test booleans; adjusted p-values stay within [0, 1]."""
    run = pipeline_runs["end_to_end"]
    table = run.assoc("body_size_mm")

    causal = table[table["variant_id"] == f"rs{CAUSAL_INDEX}"].iloc[0]
    assert bool(causal["bonferroni_sig"]) is True
    assert bool(causal["fdr_sig"]) is True
    assert causal["bonferroni_adjusted_p"] <= 1e-6
    assert table["bonferroni_sig"].sum() == (table["p_value"] <= 0.05 / len(table)).sum()
    assert (table["bonferroni_adjusted_p"] <= 1.0).all()
    assert (table["fdr_adjusted_p"] <= 1.0).all()

    report = (run.results_dir / "gwas_summary_report.txt").read_text()
    assert "Body Size (mm)" in report
    assert "HWE P-value Threshold: 1e-06" in report


def test_full_scale_summary_counts_every_trait(pipeline_runs: Dict[str, PipelineRun]) -> None:
    """The full-scale summary reports raw/QC-passing variants and the whole trait panel."""
    run = pipeline_runs["full_scale"]
    summary = json.loads((run.results_dir / "gwas_summary_full.json").read_text())

    assert summary["n_samples"] == N_SAMPLES
    assert summary["n_variants_raw"] == N_POLYMORPHIC + N_MONOMORPHIC
    assert summary["n_variants_qc_pass"] == N_POLYMORPHIC
    assert summary["traits_tested"] == len(SEVEN_TRAITS)
    assert summary["total_tests"] == N_POLYMORPHIC * len(SEVEN_TRAITS)
    by_trait = {entry["trait"]: entry for entry in summary["trait_summary"]}
    assert by_trait["Body Size"]["bonferroni_sig"] >= 1
    assert by_trait["Body Size"]["top_p_value"] < 1e-10


@pytest.mark.parametrize("name", sorted(SPECS))
def test_pipeline_exits_cleanly_when_qc_removes_every_variant(
    name: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """With nothing left after QC there is nothing to test: exit with status 1 before PCA/association."""
    config_path = build_cohort(tmp_path, n_polymorphic=0, n_monomorphic=8)

    with pytest.raises(SystemExit) as excinfo:
        SPECS[name].module.run(str(config_path))

    assert excinfo.value.code == 1
    assert "No variants passed QC" in capsys.readouterr().out
    assert not (tmp_path / "work" / SPECS[name].results_subdir / "pca_results.json").exists()


@pytest.mark.parametrize("name", sorted(SPECS))
def test_reordered_phenotype_rows_give_identical_associations(
    name: str, pipeline_runs: Dict[str, PipelineRun], tmp_path: Path
) -> None:
    """Phenotypes are matched to VCF samples by sample_id, so row order cannot change results."""
    config = build_cohort(tmp_path)
    pheno_path = tmp_path / "phenotypes.tsv"
    pheno = pd.read_csv(pheno_path, sep="\t")
    pheno.iloc[::-1].to_csv(pheno_path, sep="\t", index=False)

    spec = SPECS[name]
    spec.module.run(str(config))
    reordered = pd.read_csv(
        tmp_path / "work" / spec.results_subdir / spec.assoc_template.format(trait="body_size_mm"), sep="\t"
    )
    forward = pipeline_runs[name].assoc("body_size_mm")

    columns = ["variant_id", "beta", "p_value"]
    pd.testing.assert_frame_equal(
        reordered[columns].sort_values("variant_id").reset_index(drop=True),
        forward[columns].sort_values("variant_id").reset_index(drop=True),
    )


def _phenotype_table(ids: List[str]) -> pd.DataFrame:
    return pd.DataFrame({"sample_id": ids, "trait": [float(i) for i in range(len(ids))]})


def test_phenotypes_are_returned_in_vcf_sample_order() -> None:
    from metainformant.gwas.workflow.workflow_config import _phenotypes_in_sample_order

    assert _phenotypes_in_sample_order(_phenotype_table(["b", "a", "c"]), "trait", ["a", "b", "c"]) == [1.0, 0.0, 2.0]


@pytest.mark.parametrize(
    ("ids", "message"),
    [
        (["a", "b"], "1 missing"),
        (["a", "b", "c", "d"], "1 extra"),
        (["a", "a", "b", "c"], "duplicate sample_id"),
    ],
)
def test_phenotype_tables_that_do_not_cover_the_vcf_samples_are_refused(ids: List[str], message: str) -> None:
    from metainformant.gwas.workflow.workflow_config import _phenotypes_in_sample_order

    with pytest.raises(ValueError, match=message):
        _phenotypes_in_sample_order(_phenotype_table(ids), "trait", ["a", "b", "c"])


def test_phenotype_table_without_sample_ids_is_refused() -> None:
    from metainformant.gwas.workflow.workflow_config import _phenotypes_in_sample_order

    with pytest.raises(ValueError, match="sample_id"):
        _phenotypes_in_sample_order(pd.DataFrame({"trait": [1.0]}), "trait", ["a"])
