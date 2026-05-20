"""Across-species downstream analysis orchestrator.

Provides a unified interface for multi-species expression comparison,
divergence analysis, and ortholog expression mapping.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from metainformant.core.utils import logging
from metainformant.rna.analysis.cross_species import (
    build_ortholog_map,
    compare_expression_across_species,
    compute_expression_divergence_matrix,
    map_expression_to_orthologs,
)

logger = logging.get_logger(__name__)


class AcrossSpeciesOrchestrator:
    """Orchestrator for across-species RNA-seq comparative analysis."""

    def __init__(self, ortholog_table_path: Path, output_dir: Path):
        self.ortholog_table_path = Path(ortholog_table_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.species_expressions: Dict[str, pd.DataFrame] = {}
        self.ortholog_df: Optional[pd.DataFrame] = None

    def load_orthologs(self) -> None:
        """Load the joint ortholog table (from generate_orthologs.py)."""
        logger.info(f"Loading ortholog table from {self.ortholog_table_path}")
        self.ortholog_df = pd.read_csv(self.ortholog_table_path, sep="\t")

    def load_species_expression(self, species_name: str, abundance_path: Path) -> None:
        """Load overall expression profile for a species."""
        logger.info(f"Loading expression for {species_name} from {abundance_path}")
        expr_df = pd.read_csv(abundance_path, sep="\t", index_col=0)
        self.species_expressions[species_name] = expr_df

    def generate_pairwise_maps(self) -> Dict[Tuple[str, str], Dict[str, List[str]]]:
        """Generate all pairwise ortholog maps from the joint table."""
        if self.ortholog_df is None:
            raise ValueError("Ortholog table not loaded.")

        species_cols = [c for c in self.ortholog_df.columns if c != "Orthogroup"]
        loaded_species = list(self.species_expressions.keys())
        valid_species = [s for s in species_cols if s in loaded_species]

        maps = {}
        for sp_a in valid_species:
            for sp_b in valid_species:
                if sp_a == sp_b:
                    continue

                # Explode the comma separated lists into one-to-one or one-to-many pairs
                pairs = []
                for _, row in self.ortholog_df.iterrows():
                    a_genes = str(row[sp_a]).split(",") if pd.notna(row[sp_a]) else []
                    b_genes = str(row[sp_b]).split(",") if pd.notna(row[sp_b]) else []
                    for ag in a_genes:
                        ag = ag.strip()
                        if not ag or ag == "nan":
                            continue
                        for bg in b_genes:
                            bg = bg.strip()
                            if not bg or bg == "nan":
                                continue
                            pairs.append({sp_a: ag, sp_b: bg})

                pairs_df = pd.DataFrame(pairs)
                if not pairs_df.empty:
                    maps[(sp_a, sp_b)] = build_ortholog_map(pairs_df, sp_a, sp_b)

        return maps

    def run_comparative_analysis(self) -> None:
        """Execute full comparative suite across species."""
        if not self.species_expressions:
            raise ValueError("No species expression data loaded.")

        if self.ortholog_df is None:
            self.load_orthologs()

        logger.info("Generating pairwise ortholog mappings...")
        orth_maps = self.generate_pairwise_maps()

        logger.info("Computing global expression conservation...")
        conservation_df = compare_expression_across_species(self.species_expressions, orth_maps)

        out_cons = self.output_dir / "conservation_scores.tsv"
        conservation_df.to_csv(out_cons, sep="\t", index=False)
        logger.info(f"Saved conservation scores to {out_cons}")

        logger.info("Computing expression divergence matrix...")
        # Map everyone to the first species' ortholog space to compute divergence across a unified axis
        ref_species = list(self.species_expressions.keys())[0]
        logger.info(f"Mapping all species to reference space: {ref_species}")

        shared_expr = {ref_species: self.species_expressions[ref_species]}

        for sp in self.species_expressions:
            if sp == ref_species:
                continue
            pair_map = orth_maps.get((sp, ref_species))
            if pair_map:
                mapped = map_expression_to_orthologs(self.species_expressions[sp], pair_map)
                shared_expr[sp] = mapped

        div_matrix = compute_expression_divergence_matrix(shared_expr)
        out_div = self.output_dir / "expression_divergence_matrix.tsv"
        div_matrix.to_csv(out_div, sep="\t")
        logger.info(f"Saved expression divergence matrix to {out_div}")
