from pathlib import Path

import pandas as pd

from metainformant.epigenome import (
    read_bedgraph,
    load_cpg_table,
    compute_beta_values,
    summarize_beta_by_chromosome,
)


def test_read_bedgraph_parses_minimal_file():
    repo_root = Path(__file__).resolve().parents[1]
    bedgraph_path = repo_root / "tests/data/epigenome/example.bedgraph"
    df = read_bedgraph(bedgraph_path)
    # Expect four columns with correct names and at least 3 rows from the fixture
    assert list(df.columns) == ["chrom", "start", "end", "value"]
    assert len(df) == 3
    # Ensure dtypes are as expected
    assert pd.api.types.is_integer_dtype(df["start"]) and pd.api.types.is_integer_dtype(
        df["end"]
    )
    assert pd.api.types.is_numeric_dtype(df["value"]) and df["chrom"].dtype == object


def test_compute_methylation_beta_and_summary():
    repo_root = Path(__file__).resolve().parents[1]
    table_path = repo_root / "tests/data/epigenome/cpg_counts.tsv"
    df = load_cpg_table(table_path)
    df = compute_beta_values(df)
    # All rows should have a beta between 0 and 1 inclusive
    assert (df["beta"] >= 0).all() and (df["beta"] <= 1).all()

    summary = summarize_beta_by_chromosome(df)
    # Expect both chromosomes from the fixture to be present
    assert set(summary.index.tolist()) == {"chr1", "chr2"}
    # Means should be finite numbers
    assert pd.notnull(summary.loc["chr1"]).all()
    assert pd.notnull(summary.loc["chr2"]).all()


