"""
Expression Data Loading Module.

and handles loading of RNA-seq expression data, specifically formatted
outputs from the Amalgkit pipeline (kallisto quantification).
"""

from pathlib import Path
from typing import List, Optional, Union

import pandas as pd

from metainformant.core.utils.logging import get_logger
from metainformant.rna.core.sample_utils import find_quantification_file

logger = get_logger(__name__)

METRIC_ALIASES = {
    "tpm": ("tpm", "TPM"),
    "est_counts": ("est_counts", "NumReads", "numreads"),
    "counts": ("est_counts", "NumReads", "numreads"),
    "numreads": ("NumReads", "numreads", "est_counts"),
}
TARGET_ID_COLUMNS = ("target_id", "Name", "name")


class ExpressionLoader:
    """
    Loads and structures expression data for GWAS/eQTL analysis.
    """

    def __init__(self, work_dir: Union[str, Path]):
        """
        Initialize loader with Amalgkit working directory.

        Args:
            work_dir: Path to the main species working directory containing 'quant/'.
        """
        self.work_dir = Path(work_dir)
        self.quant_dir = self.work_dir / "quant"

    def load_amalgkit_quant(self, metric: str = "tpm", samples: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Load quantification data into a sample x transcript matrix.

        Args:
            metric: Column to extract from quantification files. Common values:
                'tpm' or 'est_counts' for kallisto, 'TPM' or 'NumReads' for salmon.
            samples: List of sample IDs to load. If None, load all found in quant_dir.

        Returns:
            pd.DataFrame: Rows = Samples, Columns = Transcripts.
        """
        if not self.quant_dir.exists():
            raise FileNotFoundError(f"Quantification directory not found: {self.quant_dir}")

        data = {}

        # Discover samples if not provided
        if samples is None:
            samples = [
                d.name
                for d in self.quant_dir.iterdir()
                if d.is_dir() and find_quantification_file(d, d.name) is not None
            ]

        logger.info(f"Loading {metric} for {len(samples)} samples...")

        for sample_id in samples:
            abundance_file = find_quantification_file(self.quant_dir / sample_id, sample_id)

            if abundance_file is None:
                logger.warning(f"Missing abundance file for {sample_id}, skipping.")
                continue

            try:
                # Load single sample
                df = pd.read_csv(abundance_file, sep="\t")

                metric_col = _resolve_column(df, metric, METRIC_ALIASES.get(metric.lower(), (metric,)))
                target_col = _resolve_column(df, "target_id", TARGET_ID_COLUMNS)

                if metric_col is None:
                    raise ValueError(f"Metric '{metric}' not found in {abundance_file}")
                if target_col is None:
                    raise ValueError(f"Transcript ID column not found in {abundance_file}")

                # Use target_id (transcript ID) as index
                # Store series in dict
                data[sample_id] = df.set_index(target_col)[metric_col]

            except Exception as e:
                logger.error(f"Failed to load {sample_id}: {e}")

        # Assemble DataFrame (Samples as rows is standard for GWAS covariates/phenotypes)
        # pd.DataFrame(data) creates Columns=Samples, Rows=Transcripts
        # Transpose to get Samples=Rows
        if not data:
            logger.warning("No data loaded.")
            return pd.DataFrame()

        # create DataFrame where columns are samples, then transpose
        files_df = pd.DataFrame(data).T
        files_df.index.name = "sample_id"

        logger.info(f"Loaded expression matrix: {files_df.shape}")
        return files_df


def _resolve_column(df: pd.DataFrame, requested: str, candidates: tuple[str, ...]) -> str | None:
    """Resolve a column name using exact or case-insensitive aliases."""
    if requested in df.columns:
        return requested
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    lower_to_column = {str(col).lower(): col for col in df.columns}
    for candidate in candidates:
        column = lower_to_column.get(candidate.lower())
        if column is not None:
            return str(column)
    return None
