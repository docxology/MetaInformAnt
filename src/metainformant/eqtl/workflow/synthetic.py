"""Synthetic and real eQTL input-data construction.

Provides generators for synthetic expression/genotype matrices with known
cis-eQTL effects (for demonstrations and zero-mocks tests) and loaders that
prepare real Amalgkit quantification data for eQTL scanning.

Example:
    >>> from metainformant.eqtl.workflow.synthetic import create_synthetic_data
    >>> expr, geno, gene_pos, var_pos = create_synthetic_data(
    ...     n_genes=5, n_variants=25, n_samples=12
    ... )
    >>> sorted(expr.shape)
    [12, 5]
"""

from __future__ import annotations

import logging
import re

import numpy as np
import pandas as pd

from metainformant.rna.core.sample_utils import find_quantification_file

logger = logging.getLogger(__name__)


def create_synthetic_data(
    n_genes: int = 100,
    n_variants: int = 500,
    n_samples: int = 50,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Create synthetic expression and genotype data with known eQTL structure.

    Simulates realistic eQTL patterns including cis-eQTLs with strong effects
    and background noise variants. Roughly 30% of genes receive a true
    cis-eQTL effect on their first cis-window variant.

    Args:
        n_genes: Number of genes (must satisfy n_variants >= 5 * n_genes).
        n_variants: Total number of variants (5 per gene in cis window).
        n_samples: Number of samples.
        seed: Random seed for reproducibility.

    Returns:
        Tuple of (expression_matrix, genotype_matrix, gene_positions,
        variant_positions). Expression is genes x samples; genotypes are
        variants x samples with dosages in {0, 1, 2}.
    """
    min_variants = 5 * n_genes
    if n_variants < min_variants:
        raise ValueError(f"n_variants ({n_variants}) must be >= 5 * n_genes ({min_variants})")
    if n_samples < 2:
        raise ValueError(f"n_samples ({n_samples}) must be >= 2")

    rng = np.random.default_rng(seed)
    logger.info(f"Generating synthetic data: {n_genes} genes, {n_variants} variants, " f"{n_samples} samples")

    # Sample IDs
    sample_ids = [f"SRR{10000000 + i}" for i in range(n_samples)]

    # Gene positions (spread across chromosomes 1-5)
    gene_ids = [f"LOC{100000 + i}" for i in range(n_genes)]
    chromosomes = [str((i % 5) + 1) for i in range(n_genes)]
    tss_positions = [1_000_000 + (i * 50_000) for i in range(n_genes)]

    gene_positions = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "chrom": chromosomes,
            "tss_position": tss_positions,
        }
    )

    # Variant positions (5 per gene in cis window)
    variant_ids = []
    var_chromosomes = []
    var_positions = []

    for gene_idx in range(n_genes):
        for var_offset in range(5):
            variant_ids.append(f"rs{gene_idx * 5 + var_offset}")
            var_chromosomes.append(chromosomes[gene_idx])
            # Position within cis window
            var_positions.append(tss_positions[gene_idx] + (var_offset - 2) * 10_000)

    variant_positions = pd.DataFrame(
        {
            "variant_id": variant_ids,
            "chrom": var_chromosomes,
            "position": var_positions,
        }
    )

    # Genotypes: dosages 0, 1, 2 with MAF ~ 0.3
    genotypes = rng.choice([0, 1, 2], size=(len(variant_ids), n_samples), p=[0.5, 0.35, 0.15])
    genotype_matrix = pd.DataFrame(genotypes, index=variant_ids, columns=sample_ids)

    # Expression: baseline + eQTL effect for ~30% of genes
    expression = np.zeros((n_genes, n_samples))

    for gene_idx in range(n_genes):
        # Baseline expression (log-scale, ~5-10)
        baseline = rng.uniform(5, 10)
        noise = rng.normal(0, 0.5, n_samples)

        # Add eQTL effect for first variant of each gene (30% chance)
        if rng.random() < 0.3:
            var_idx = gene_idx * 5  # First variant for this gene
            eqtl_effect = rng.uniform(0.5, 2.0)
            effect = genotypes[var_idx] * eqtl_effect
        else:
            effect = 0

        expression[gene_idx] = baseline + effect + noise

    expression_matrix = pd.DataFrame(expression, index=gene_ids, columns=sample_ids)

    logger.info(f"Created expression matrix: {expression_matrix.shape}")
    logger.info(f"Created genotype matrix: {genotype_matrix.shape}")

    return expression_matrix, genotype_matrix, gene_positions, variant_positions


def parse_gene_positions(gene_ids: list[str]) -> pd.DataFrame:
    """Parse approximate gene positions from kallisto target IDs.

    Target format: ``lcl|NC_037638.1_mrna_XM_623972.6_1``. The chromosome is
    extracted from the NC_* accession prefix; positions are synthetic (spread
    evenly) because real coordinates require a GFF annotation.

    Args:
        gene_ids: Transcript/gene target identifiers.

    Returns:
        DataFrame with columns [gene_id, chrom, tss_position].
    """
    positions = []

    # Chromosome accession pattern: NC_037638.1 etc., possibly prefixed by lcl|
    chrom_pattern = re.compile(r"(NC_\d+\.\d+|NW_\d+\.\d+|NT_\d+\.\d+)")

    for idx, gid in enumerate(gene_ids):
        match = chrom_pattern.search(gid.replace("lcl|", ""))
        chrom = match.group(1) if match else "unknown"

        # Spread across genome (~250Mb / n_genes); real positions need a GFF
        position = 1_000_000 + (idx * 10_000)

        positions.append(
            {
                "gene_id": gid,
                "chrom": chrom,
                "tss_position": position,
            }
        )

    return pd.DataFrame(positions)


def create_synthetic_genotypes(
    sample_ids: list[str],
    gene_positions: pd.DataFrame,
    variants_per_gene: int = 1,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create synthetic genotype data matched to gene positions.

    Note: in real analysis genotypes come from matched WGS/genotyping data;
    this generator exists for demonstration and workflow testing.

    Args:
        sample_ids: Sample identifiers (genotype columns).
        gene_positions: DataFrame with columns [gene_id, chrom, tss_position].
        variants_per_gene: Number of cis variants simulated per gene.
        seed: Random seed for reproducibility.

    Returns:
        Tuple of (genotype_matrix, variant_positions). Genotypes are
        variants x samples dosages in {0, 1, 2} with MAF ~ 0.25.
    """
    logger.info("Creating synthetic genotypes (real genotypes not available)")
    rng = np.random.default_rng(seed)

    variant_ids = []
    var_chroms = []
    var_positions = []

    for _, row in gene_positions.iterrows():
        for i in range(variants_per_gene):
            var_id = f"var_{row['gene_id'][:20]}_{i}"
            variant_ids.append(var_id)
            var_chroms.append(row["chrom"])
            var_positions.append(int(row["tss_position"]) + (i - 1) * 5000)

    variant_positions = pd.DataFrame(
        {
            "variant_id": variant_ids,
            "chrom": var_chroms,
            "position": var_positions,
        }
    )

    # Generate dosages with MAF ~ 0.25
    n_variants = len(variant_ids)
    n_samples = len(sample_ids)
    genotypes = rng.choice([0, 1, 2], size=(n_variants, n_samples), p=[0.56, 0.32, 0.12])

    genotype_matrix = pd.DataFrame(genotypes, index=variant_ids, columns=sample_ids)

    return genotype_matrix, variant_positions


def load_real_expression_data(
    quant_dir,
    max_samples: int = 100,
    min_tpm: float = 1.0,
    max_genes: int = 1000,
) -> tuple[pd.DataFrame, list[str]]:
    """Load real expression data from recognized per-sample quant files.

    Args:
        quant_dir: Directory containing sample subdirectories (Path-like).
        max_samples: Maximum samples to load (for speed).
        min_tpm: Minimum mean TPM to include a gene.
        max_genes: Maximum genes to include (top by mean expression).

    Returns:
        Expression matrix (genes x samples) and list of sample IDs.

    Raises:
        ValueError: If no expression data could be loaded from quant_dir.
    """
    from pathlib import Path

    quant_path = Path(quant_dir)
    logger.info(f"Loading real expression data from {quant_path}")

    sample_dirs = sorted([d for d in quant_path.iterdir() if d.is_dir()])
    logger.info(f"Found {len(sample_dirs)} sample directories")

    # Limit samples for speed
    sample_dirs = sample_dirs[:max_samples]
    logger.info(f"Loading {len(sample_dirs)} samples")

    expression_data: dict[str, pd.Series] = {}
    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        abundance_file = find_quantification_file(sample_dir, sample_id)
        if abundance_file is None:
            continue

        df = pd.read_csv(abundance_file, sep="\t")
        target_col = "target_id" if "target_id" in df.columns else "Name"
        tpm_col = "tpm" if "tpm" in df.columns else "TPM"
        if target_col not in df.columns or tpm_col not in df.columns:
            logger.warning(f"Skipping {sample_id}: missing transcript or TPM column in {abundance_file}")
            continue

        expression_data[sample_id] = df.set_index(target_col)[tpm_col]

    if not expression_data:
        raise ValueError("No expression data loaded")

    # Build expression matrix
    expr_matrix = pd.DataFrame(expression_data)
    logger.info(f"Loaded expression matrix: {expr_matrix.shape}")

    # Filter low-expression genes
    mean_tpm = expr_matrix.mean(axis=1)
    expressed_genes = mean_tpm[mean_tpm >= min_tpm].index
    expr_matrix = expr_matrix.loc[expressed_genes]
    logger.info(f"After filtering (TPM >= {min_tpm}): {len(expr_matrix)} genes")

    # Keep only top genes by expression for speed
    if len(expr_matrix) > max_genes:
        top_genes = mean_tpm.loc[expr_matrix.index].nlargest(max_genes).index
        expr_matrix = expr_matrix.loc[top_genes]
        logger.info(f"Selected top {max_genes} genes by expression")

    return expr_matrix, list(expression_data.keys())
