from __future__ import annotations

from pathlib import Path

import pandas as pd

from metainformant.core.io import read_delimited


def load_cpg_table(path: str | Path) -> pd.DataFrame:
    """Load a simple CpG count table with columns chrom,pos,methylated,unmethylated (TSV).

    Returns a DataFrame with proper dtypes: chrom as object, pos as int, counts as int.
    """
    rows = list(read_delimited(path, delimiter="\t"))
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=["chrom", "pos", "methylated", "unmethylated"])  # type: ignore[return-value]
    df["pos"] = df["pos"].astype(int)
    df["methylated"] = df["methylated"].astype(int)
    df["unmethylated"] = df["unmethylated"].astype(int)
    return df


def compute_beta_values(df: pd.DataFrame) -> pd.DataFrame:
    """Compute beta = methylated / (methylated + unmethylated), safe for zeros.

    Returns a new DataFrame with an added 'beta' column.
    """
    total = (df["methylated"].astype(float) + df["unmethylated"].astype(float))
    # Avoid division by zero; where total is 0, define beta as 0.0
    beta = pd.Series(0.0, index=df.index)
    nonzero = total > 0
    beta.loc[nonzero] = (df.loc[nonzero, "methylated"].astype(float) / total.loc[nonzero])
    result = df.copy()
    result["beta"] = beta
    return result


def summarize_beta_by_chromosome(df_with_beta: pd.DataFrame) -> pd.DataFrame:
    """Return per-chromosome mean beta as a DataFrame indexed by chrom.

    A single column 'beta_mean' is returned so that downstream code can
    select a row (e.g., summary.loc["chr1"]) yielding a Series, which
    plays well with vectorized assertions in tests.
    Requires a 'beta' column.
    """
    if "beta" not in df_with_beta.columns:
        raise ValueError("beta column missing; compute_beta_values first")
    means = df_with_beta.groupby("chrom", as_index=True)["beta"].mean()
    return means.to_frame(name="beta_mean")



