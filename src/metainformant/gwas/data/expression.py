"""
Expression Data Loading Module.

and handles loading of RNA-seq expression data, specifically formatted
outputs from the Amalgkit pipeline (kallisto quantification).
"""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


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

    def load_amalgkit_quant(
        self, 
        metric: str = "tpm", 
        samples: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Load quantification data into a sample x transcript matrix.

        Args:
            metric: Column to extract from abundance.tsv ('tpm' or 'est_counts').
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
                d.name for d in self.quant_dir.iterdir() 
                if d.is_dir() and (d / "abundance.tsv").exists()
            ]
        
        logger.info(f"Loading {metric} for {len(samples)} samples...")
        
        for sample_id in samples:
            abundance_file = self.quant_dir / sample_id / "abundance.tsv"
            
            if not abundance_file.exists():
                logger.warning(f"Missing abundance file for {sample_id}, skipping.")
                continue
                
            try:
                # Load single sample
                df = pd.read_csv(abundance_file, sep="\t")
                
                if metric not in df.columns:
                    raise ValueError(f"Metric '{metric}' not found in {abundance_file}")
                
                # Use target_id (transcript ID) as index
                # Store series in dict
                data[sample_id] = df.set_index("target_id")[metric]
                
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
