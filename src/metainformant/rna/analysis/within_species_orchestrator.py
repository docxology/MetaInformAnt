"""Within-species downstream analysis orchestrator.

Provides a unified interface for within-species differential expression
and principal component analysis, processing amalgkit matrix outputs.
"""

from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import pandas as pd
from pathlib import Path

from metainformant.core.utils import logging
from metainformant.rna.analysis.expression_analysis import (
    differential_expression,
    pca_analysis,
    prepare_volcano_data
)

logger = logging.get_logger(__name__)

class WithinSpeciesOrchestrator:
    """Orchestrator for within-species RNA-seq analysis."""

    def __init__(self, species_name: str, abundance_path: Path, metadata_path: Path, output_dir: Path):
        self.species_name = species_name
        self.abundance_path = Path(abundance_path)
        self.metadata_path = Path(metadata_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.counts_df: Optional[pd.DataFrame] = None
        self.metadata_df: Optional[pd.DataFrame] = None

    def load_data(self) -> None:
        """Load abundance and metadata matrices."""
        logger.info(f"[{self.species_name}] Loading abundance matrix: {self.abundance_path}")
        self.counts_df = pd.read_csv(self.abundance_path, sep="\t", index_col=0)
        
        logger.info(f"[{self.species_name}] Loading metadata: {self.metadata_path}")
        self.metadata_df = pd.read_csv(self.metadata_path, sep="\t")
        
        # Ensure sample IDs match (amalgkit might use 'run' or 'sample')
        sample_col = "run" if "run" in self.metadata_df.columns else self.metadata_df.columns[0]
        self.metadata_df.set_index(sample_col, inplace=True)
        
        # Align columns
        common_samples = list(set(self.counts_df.columns) & set(self.metadata_df.index))
        if not common_samples:
            raise ValueError("No common samples between abundance matrix and metadata.")
            
        self.counts_df = self.counts_df[common_samples]
        self.metadata_df = self.metadata_df.loc[common_samples]
        logger.info(f"[{self.species_name}] Aligned {len(common_samples)} samples.")

    def run_pca(self, n_components: int = 2) -> Dict[str, pd.DataFrame]:
        """Run PCA on the expression data."""
        if self.counts_df is None:
            raise ValueError("Data not loaded. Call load_data() first.")
            
        logger.info(f"[{self.species_name}] Running PCA ({n_components} components)...")
        # Log2 transform counts for basic PCA scaling
        log_counts = np.log2(self.counts_df + 1)
        pca_res = pca_analysis(log_counts, n_components=n_components)
        
        out_path = self.output_dir / f"{self.species_name}_pca_coordinates.tsv"
        pca_res["transformed"].to_csv(out_path, sep="\t")
        logger.info(f"[{self.species_name}] Saved PCA coordinates to {out_path}")
        
        return pca_res

    def run_differential_expression(self, condition_col: str, reference: Optional[str] = None) -> Optional[pd.DataFrame]:
        """Run standard pairwise DE comparing exactly two conditions in a metadata column."""
        if self.counts_df is None or self.metadata_df is None:
            raise ValueError("Data not loaded. Call load_data() first.")
            
        if condition_col not in self.metadata_df.columns:
            logger.warning(f"[{self.species_name}] Condition '{condition_col}' not found in metadata. Skipping DE.")
            return None
            
        # Filter samples that have a valid value for this condition
        valid_meta = self.metadata_df[self.metadata_df[condition_col].notna()]
        unique_vals = valid_meta[condition_col].unique()
        
        if len(unique_vals) != 2:
            logger.warning(f"[{self.species_name}] DE requires exactly 2 conditions for '{condition_col}'. Found {len(unique_vals)}: {unique_vals}. Skipping.")
            return None
            
        valid_samples = valid_meta.index.tolist()
        subset_counts = self.counts_df[valid_samples]
        conditions = valid_meta[condition_col]
        
        logger.info(f"[{self.species_name}] Running DE on {condition_col}: {unique_vals}")
        de_res = differential_expression(
            subset_counts, 
            conditions, 
            method="ttest",
            reference=reference
        )
        
        volcano_data = prepare_volcano_data(de_res)
        out_path = self.output_dir / f"{self.species_name}_DE_{condition_col}.tsv"
        volcano_data.to_csv(out_path, sep="\t", index=False)
        logger.info(f"[{self.species_name}] Saved DE results to {out_path}")
        
        return volcano_data

    def run_all(self, condition_cols: List[str] = ["tissue", "sex", "caste", "developmental_stage"]) -> None:
        """Execute full within-species analysis suite."""
        self.load_data()
        self.run_pca()
        
        for col in condition_cols:
            try:
                self.run_differential_expression(col)
            except Exception as e:
                logger.error(f"[{self.species_name}] DE failed on {col}: {e}")
