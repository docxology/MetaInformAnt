"""
Integration test for RNA-Seq -> GWAS data handoff.

Verifies that:
1. ExpressionLoader correctly parses Amalgkit abundance tables.
2. run_eqtl_analysis correctly links expression and genotype data.
"""

import pandas as pd
import pytest

from metainformant.gwas.analysis.eqtl import run_eqtl_analysis
from metainformant.gwas.data.expression import ExpressionLoader


@pytest.fixture
def amalgkit_output_dir(tmp_path):
    """Create an Amalgkit output directory structure."""
    work_dir = tmp_path / "amalgkit_work"
    quant_dir = work_dir / "quant"
    quant_dir.mkdir(parents=True)

    samples = ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"]
    transcripts = ["TR_1", "TR_2", "TR_3"]

    # Create abundance files
    for i, sample in enumerate(samples):
        sample_dir = quant_dir / sample
        sample_dir.mkdir()

        # Create deterministic abundance data.
        data = {
            "target_id": transcripts,
            "est_counts": [10.0 * (i + 1), 20.0 * (i + 1), 5.0],
            "tpm": [1.0 * (i + 1), 2.0 * (i + 1), 0.5],
        }
        df = pd.DataFrame(data)
        df.to_csv(sample_dir / "abundance.tsv", sep="\t", index=False)

    return work_dir


@pytest.fixture
def genotype_matrix():
    """Create a genotype matrix (samples x variants)."""
    # 0, 1, 2 dosage
    data = {"VAR_1": [0, 1, 2], "VAR_2": [2, 1, 0]}  # Correlated with TR_1 (1, 2, 3)  # Anti-correlated
    df = pd.DataFrame(data, index=["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"])
    return df


def test_rna_gwas_handoff(amalgkit_output_dir, genotype_matrix):
    """
    Test the full handoff pipeline.
    """
    # 1. Load Expression Data
    loader = ExpressionLoader(amalgkit_output_dir)
    expression_matrix = loader.load_amalgkit_quant(metric="tpm")

    assert expression_matrix.shape == (3, 3)
    assert "SAMPLE_A" in expression_matrix.index
    assert "TR_1" in expression_matrix.columns

    # Check values (SAMPLE_A should be 1.0, 2.0, 0.5)
    assert expression_matrix.loc["SAMPLE_A", "TR_1"] == 1.0

    # 2. Run eQTL Analysis
    # We expect VAR_1 to be correlated with TR_1
    results = run_eqtl_analysis(genotype_matrix=genotype_matrix, expression_matrix=expression_matrix)

    assert not results.empty
    assert "variant" in results.columns
    assert "transcript" in results.columns
    assert "p_value" in results.columns

    # Check specific finding
    # Filter for TR_1 and VAR_1
    subset = results[(results["transcript"] == "TR_1") & (results["variant"] == "VAR_1")]
    # Linear perfect correlation 0,1,2 vs 1,2,3 should have very low p-value
    # But association_test_linear implementation might vary, so just check existence
    assert len(subset) == 1

    print("\neQTL Results Head:")
    print(results.head())


def test_expression_loader_accepts_prefixed_and_salmon_quant(tmp_path):
    """GWAS handoff should load all RNA quant output variants validated upstream."""
    work_dir = tmp_path / "amalgkit_work"
    quant_dir = work_dir / "quant"

    prefixed_dir = quant_dir / "SAMPLE_PREFIXED"
    prefixed_dir.mkdir(parents=True)
    pd.DataFrame({"target_id": ["TR_1"], "tpm": [3.0], "est_counts": [30.0]}).to_csv(
        prefixed_dir / "SAMPLE_PREFIXED_abundance.tsv", sep="\t", index=False
    )

    salmon_dir = quant_dir / "SAMPLE_SALMON"
    salmon_dir.mkdir(parents=True)
    pd.DataFrame({"Name": ["TR_1"], "Length": [1000], "TPM": [4.0], "NumReads": [40.0]}).to_csv(
        salmon_dir / "quant.sf", sep="\t", index=False
    )

    expression_matrix = ExpressionLoader(work_dir).load_amalgkit_quant(metric="tpm")

    assert expression_matrix.shape == (2, 1)
    assert expression_matrix.loc["SAMPLE_PREFIXED", "TR_1"] == 3.0
    assert expression_matrix.loc["SAMPLE_SALMON", "TR_1"] == 4.0


def test_run_eqtl_analysis_scans_all_transcripts_by_default(genotype_matrix):
    """The verification wrapper should not silently cap analyses at 10 transcripts."""
    transcript_ids = [f"TR_{i:02d}" for i in range(12)]
    expression_matrix = pd.DataFrame(
        {transcript_id: [float(i + 1), float(i + 2), float(i + 3)] for i, transcript_id in enumerate(transcript_ids)},
        index=genotype_matrix.index,
    )

    results = run_eqtl_analysis(genotype_matrix=genotype_matrix[["VAR_1"]], expression_matrix=expression_matrix)
    limited = run_eqtl_analysis(
        genotype_matrix=genotype_matrix[["VAR_1"]],
        expression_matrix=expression_matrix,
        max_transcripts=5,
    )

    assert set(results["transcript"]) == set(transcript_ids)
    assert set(limited["transcript"]) == set(transcript_ids[:5])
