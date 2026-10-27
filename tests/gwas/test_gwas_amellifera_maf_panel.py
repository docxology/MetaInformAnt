"""Regression tests: the A. mellifera pipeline plots real minor allele frequencies, never |beta|.

Defect: ``amellifera_pipeline.generate_visualizations`` built the "MAF distribution" panel from
``min(abs(beta), 0.5)`` of every result, with a comment claiming "real MAFs from genotype data". |beta| is an effect
size, so the panel (titled "Minor Allele Frequency Distribution", with "Rare variants" and "Common variants" counts)
described a quantity that has nothing to do with allele frequencies. The association results did not carry a MAF at
all (only the summary tables did, computed by the writers), so the panel could not have been drawn honestly.

Now ``run_gwas`` and ``execute_gwas_workflow`` attach to every result the minor allele frequency of the genotypes its
test ran on (``maf``), and the panel is drawn from those; results without one are left out, and when none has one the
panel is skipped with a warning.

Oracles: allele counts worked out by hand from the dosages; rendering is deterministic, so panels that must show the
same data are byte-identical and panels that must differ are not.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, List

import pytest

from metainformant.gwas.workflow.amellifera_pipeline import generate_visualizations
from metainformant.gwas.workflow.workflow_execution import run_gwas

N = 12


def results(maf_scale: float = 1.0, beta_scale: float = 1.0, with_maf: bool = True) -> List[Dict[str, Any]]:
    """Twelve results whose |beta| (0.40 to 0.51) and MAF (0.05 to 0.38 times ``maf_scale``) are far apart."""
    out = []
    for i in range(N):
        row: Dict[str, Any] = {
            "variant_id": f"v{i}",
            "chrom": "NC_037638.1" if i < 7 else "NC_037639.1",
            "pos": 1000 * (i + 1),
            "beta": beta_scale * (0.40 + 0.01 * i),
            "se": 0.1,
            "p_value": 10.0 ** -(1 + i % 6),
        }
        if with_maf:
            row["maf"] = maf_scale * (0.05 + 0.03 * i)
        out.append(row)
    return out


@pytest.fixture(scope="module", autouse=True)
def settled_plot_style(tmp_path_factory: pytest.TempPathFactory) -> None:
    """The first plots of a process set matplotlib's style globally; draw once so later runs are comparable."""
    generate_visualizations({"results": {"association_results": results()}}, tmp_path_factory.mktemp("warm_up"))


def run(rows: List[Dict[str, Any]], directory: Path) -> Dict[str, str]:
    return generate_visualizations({"results": {"association_results": rows}}, directory)


def panel(directory: Path) -> bytes:
    return (directory / "plots" / "maf_distribution.png").read_bytes()


# ---------------------------------------------------------------------------
# generate_visualizations
# ---------------------------------------------------------------------------


def test_the_panel_follows_the_mafs_of_the_results_and_not_their_effect_sizes(tmp_path: Path) -> None:
    baseline = run(results(), tmp_path / "baseline")
    run(results(beta_scale=0.3), tmp_path / "other_betas")  # the same MAFs, effect sizes a third as large
    run(results(maf_scale=0.5), tmp_path / "other_mafs")  # the same effect sizes, half the allele frequencies

    assert "maf_distribution" in baseline
    assert panel(tmp_path / "baseline") == panel(tmp_path / "other_betas")  # |beta| is not what is drawn
    assert panel(tmp_path / "baseline") != panel(tmp_path / "other_mafs")  # the MAFs are


def test_results_without_a_maf_get_no_panel_and_a_warning(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.WARNING):
        written = run(results(with_maf=False), tmp_path)

    assert "maf_distribution" not in written
    assert not (tmp_path / "plots" / "maf_distribution.png").exists()
    assert "MAF distribution skipped: the association results carry no minor allele frequency" in caplog.text
    assert {"manhattan", "qq", "volcano"} <= set(written)  # the other panels are unaffected


def test_results_with_a_maf_for_only_some_variants_plot_those_and_say_so(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    partial = results()
    for row in partial[8:]:
        del row["maf"]

    with caplog.at_level(logging.WARNING):
        run(partial, tmp_path / "partial")
    run(partial[:8], tmp_path / "first_eight")  # the eight that carry a MAF, on their own

    assert "MAF distribution drawn from 8 of 12 results" in caplog.text
    assert panel(tmp_path / "partial") == panel(tmp_path / "first_eight")


def test_a_non_finite_maf_counts_as_missing(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    rows = results()
    rows[0]["maf"] = float("nan")
    rows[1]["maf"] = None

    with caplog.at_level(logging.WARNING):
        run(rows, tmp_path / "gaps")
    run(rows[2:], tmp_path / "rest")

    assert "MAF distribution drawn from 10 of 12 results" in caplog.text
    assert panel(tmp_path / "gaps") == panel(tmp_path / "rest")


# ---------------------------------------------------------------------------
# The workflow results carry the MAF
# ---------------------------------------------------------------------------

GENOTYPES = {
    "rsA": [0, 1, 2, 0, 1, 2, 0, 1, 2, 1, 0, 2],  # 12 alt alleles of 24: MAF 0.5
    "rsB": [2, 2, 1, 0, 0, 1, 2, 1, 0, 0, 1, 1],  # 13 of 24: alt 0.5417, MAF 0.4583
    "rsC": [0, 0, 1, 1, 2, 2, 1, 0, 2, 1, 2, 0],  # 12 of 24: MAF 0.5
    "rsM": [0] * 12,  # monomorphic: MAF 0, and its test fails
    "rsD": [1, 0, 2, 1, 0, 1, 2, 2, 0, 1, 1, 0],  # 11 of 24: MAF 0.4583
}
PHENOTYPE = [1.2, 2.3, 2.9, 1.8, 3.6, 2.1, 3.9, 1.1, 2.6, 3.1, 1.9, 2.4]


def write_dataset(directory: Path, n_phenotyped: int = 12) -> Dict[str, Path]:
    calls = {0: "0/0", 1: "0/1", 2: "1/1"}
    samples = [f"S{i + 1}" for i in range(N)]
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    for i, (snp, row) in enumerate(GENOTYPES.items()):
        fields = ["chr1", str(100 * (i + 1)), snp, "A", "G", "60", "PASS", ".", "GT"]
        lines.append("\t".join(fields + [calls[g] for g in row]))
    vcf = directory / "data.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    phenotypes = directory / "phenotypes.tsv"
    body = [f"S{i + 1}\t{y}" for i, y in enumerate(PHENOTYPE[:n_phenotyped])]
    phenotypes.write_text("\n".join(["sample_id\tyield"] + body) + "\n")
    return {"vcf": vcf, "phenotypes": phenotypes}


def counted_maf(genotypes: List[int]) -> float:
    alt = sum(genotypes) / (2 * len(genotypes))
    return min(alt, 1 - alt)


def run_workflow(files: Dict[str, Path], out: Path) -> List[Dict[str, Any]]:
    config = {
        "model": "linear",
        "trait": "yield",
        "quality_control": {"min_maf": 0.0, "max_missing": 1.0, "min_hwe_p": 0.0},
        "significance_threshold": 1.0,
    }
    result = run_gwas(files["vcf"], files["phenotypes"], config, output_dir=out)
    assert result["status"] == "success"
    return result["results"]["association_results"]


def test_run_gwas_results_carry_the_maf_of_the_genotypes_their_tests_used(tmp_path: Path) -> None:
    rows = run_workflow(write_dataset(tmp_path), tmp_path / "out")

    assert [row["maf"] for row in rows] == pytest.approx([counted_maf(g) for g in GENOTYPES.values()])
    assert rows[3]["maf"] == 0.0 and rows[3]["status"] == "error"  # a failed test still has an allele frequency


def test_run_gwas_counts_only_the_samples_its_tests_used(tmp_path: Path) -> None:
    """Phenotypes for eight of the twelve samples: the tests, and so the MAFs, use the first eight."""
    rows = run_workflow(write_dataset(tmp_path, n_phenotyped=8), tmp_path / "out")

    expected = [counted_maf(g[:8]) for g in GENOTYPES.values()]
    assert expected != pytest.approx([counted_maf(g) for g in GENOTYPES.values()])  # the two sample sets differ
    assert [row["maf"] for row in rows] == pytest.approx(expected)


def test_the_written_table_and_the_rows_agree(tmp_path: Path) -> None:
    rows = run_workflow(write_dataset(tmp_path), tmp_path / "out")

    table = (tmp_path / "out" / "summary_statistics.tsv").read_text().strip().split("\n")
    header = table[0].split("\t")
    written = [float(dict(zip(header, line.split("\t")))["MAF"]) for line in table[1:]]
    assert written == pytest.approx([row["maf"] for row in rows], abs=5.1e-5)  # the table keeps four decimals


# ---------------------------------------------------------------------------
# Attaching the MAF
# ---------------------------------------------------------------------------


def test_an_existing_maf_is_kept_and_a_missing_one_is_not_invented() -> None:
    from metainformant.gwas.workflow.workflow_execution import _attach_maf

    rows: List[Dict[str, Any]] = [{"maf": 0.123}, {"MAF": 0.2}, {"maf": None}, {"maf": float("nan")}, {}]
    genotypes = [[0, 1, 2, 0], [0, 1, 2, 0], [1, 1, 0, 0], [-1, -1, -1, -1], [0, 0, 1, 1]]

    _attach_maf(rows, genotypes)

    assert rows[0] == {"maf": 0.123}  # kept
    assert rows[1] == {"MAF": 0.2}  # an upper-case MAF is kept as it is
    assert rows[2]["maf"] == pytest.approx(0.25)  # None is replaced by the counted frequency (2 of 8)
    assert "maf" not in rows[3] or math.isnan(rows[3]["maf"])  # no called genotype: nothing is invented
    assert rows[4]["maf"] == pytest.approx(0.25)


def test_rows_and_genotypes_that_do_not_line_up_get_nothing_and_a_warning(caplog: pytest.LogCaptureFixture) -> None:
    from metainformant.gwas.workflow.workflow_execution import _attach_maf

    rows: List[Dict[str, Any]] = [{}, {}, {}]

    with caplog.at_level(logging.WARNING):
        _attach_maf(rows, [[0, 1, 2], [1, 1, 1]])

    assert rows == [{}, {}, {}]
    assert "3 results for 2 genotype rows" in caplog.text


@pytest.mark.parametrize(
    ("dosages", "expected"),
    [
        ([0, 1, 2, 0, 1, 2], 0.5),  # 6 alternate alleles of 12
        ([0, 0, 0, 1], 0.125),  # 1 of 8
        ([2, 2, 2, 1], 0.125),  # 7 of 8: the alternate allele is the major one, so the MAF is 1 - 7/8
        ([0, 1, -1, float("nan"), 2], 0.5),  # only the three called dosages count: 3 alternate alleles of 6
    ],
)
def test_the_minor_allele_frequency_counts_called_dosages_only(dosages: List[float], expected: float) -> None:
    from metainformant.gwas.analysis.summary_stats import minor_allele_frequency

    assert minor_allele_frequency(dosages) == pytest.approx(expected)


def test_a_row_without_a_called_genotype_has_no_maf() -> None:
    from metainformant.gwas.analysis.summary_stats import minor_allele_frequency

    assert minor_allele_frequency([-1, -1, float("nan")]) is None
    assert minor_allele_frequency([]) is None
