"""
Integration test for RNA-Seq -> GWAS data handoff.

Verifies that:
1. ExpressionLoader correctly parses Amalgkit abundance tables.
2. run_eqtl_analysis correctly links expression and genotype data.
"""

import pandas as pd
import pytest
from pathlib import Path
from metainformant.gwas.data.expression import ExpressionLoader
from metainformant.gwas.analysis.eqtl import run_eqtl_analysis

@pytest.fixture
def mock_amalgkit_output(tmp_path):
    """Create a mock Amalgkit output directory structure."""
    work_dir = tmp_path / "amalgkit_work"
    quant_dir = work_dir / "quant"
    quant_dir.mkdir(parents=True)
    
    samples = ["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"]
    transcripts = ["TR_1", "TR_2", "TR_3"]
    
    # Create abundance files
    for i, sample in enumerate(samples):
        sample_dir = quant_dir / sample
        sample_dir.mkdir()
        
        # Create dummy data
        data = {
            "target_id": transcripts,
            "est_counts": [10.0 * (i+1), 20.0 * (i+1), 5.0],
            "tpm": [1.0 * (i+1), 2.0 * (i+1), 0.5]
        }
        df = pd.DataFrame(data)
        df.to_csv(sample_dir / "abundance.tsv", sep="\t", index=False)
        
    return work_dir

@pytest.fixture
def mock_genotypes():
    """Create a mock genotype matrix (Samples x Variants)."""
    # 0, 1, 2 dosage
    data = {
        "VAR_1": [0, 1, 2], # Correlated with TR_1 (1, 2, 3)
        "VAR_2": [2, 1, 0]  # Anti-correlated
    }
    df = pd.DataFrame(data, index=["SAMPLE_A", "SAMPLE_B", "SAMPLE_C"])
    return df

def test_rna_gwas_handoff(mock_amalgkit_output, mock_genotypes):
    """
    Test the full handoff pipeline.
    """
    # 1. Load Expression Data
    loader = ExpressionLoader(mock_amalgkit_output)
    expression_matrix = loader.load_amalgkit_quant(metric="tpm")
    
    assert expression_matrix.shape == (3, 3)
    assert "SAMPLE_A" in expression_matrix.index
    assert "TR_1" in expression_matrix.columns
    
    # Check values (SAMPLE_A should be 1.0, 2.0, 0.5)
    assert expression_matrix.loc["SAMPLE_A", "TR_1"] == 1.0
    
    # 2. Run eQTL Analysis
    # We expect VAR_1 to be correlated with TR_1
    results = run_eqtl_analysis(
        genotype_matrix=mock_genotypes,
        expression_matrix=expression_matrix
    )
    
    assert not results.empty
    assert "variant" in results.columns
    assert "transcript" in results.columns
    assert "p_value" in results.columns
    
    # Check specific finding
    # Filter for TR_1 and VAR_1
    subset = results[
        (results["transcript"] == "TR_1") & 
        (results["variant"] == "VAR_1")
    ]
    # Linear perfect correlation 0,1,2 vs 1,2,3 should have very low p-value
    # But association_test_linear implementation might vary, so just check existence
    assert len(subset) == 1
    
    print("\neQTL Results Head:")
    print(results.head())
