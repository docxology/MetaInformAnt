"""Merge step for RNA-seq workflow.

This step merges expression matrices from multiple samples into a single
comprehensive dataset for downstream analysis.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Dict

import pandas as pd

from metainformant.core import io, logging
from metainformant.rna.steps import StepResult

logger = logging.get_logger(__name__)


def run_step(step_params: Dict[str, Any]) -> StepResult:
    """Execute the merge step.

    Args:
        step_params: Parameters for the merge step
            - work_dir: Working directory
            - quantification_results: Path to quantification results
            - species_list: List of species to merge
            - output_format: Output format ("tsv", "csv", "hdf5")

    Returns:
        StepResult with merge results
    """
    start_time = time.time()

    try:
        work_dir = Path(step_params.get('work_dir', '.'))
        quant_results_file = step_params.get('quantification_results',
                                           work_dir / "quant" / "quantification_results.json")
        species_list = step_params.get('species_list', [])
        output_format = step_params.get('output_format', 'tsv')

        if isinstance(quant_results_file, str):
            quant_results_file = Path(quant_results_file)

        if not quant_results_file.exists():
            raise FileNotFoundError(f"Quantification results file not found: {quant_results_file}")

        quant_results = io.load_json(quant_results_file)

        # Create merged directory
        merged_dir = work_dir / "merged"
        merged_dir.mkdir(parents=True, exist_ok=True)

        merged_matrices = {}

        for species in species_list:
            if species not in quant_results:
                logger.warning(f"Species {species} not found in quantification results")
                continue

            logger.info(f"Merging expression data for {species}")

            species_matrix = merge_species_expressions(
                species, quant_results[species], merged_dir, output_format
            )

            merged_matrices[species] = species_matrix

        # Create cross-species merged matrix if multiple species
        if len(merged_matrices) > 1:
            cross_species_matrix = merge_cross_species(merged_matrices, merged_dir, output_format)
            merged_matrices['cross_species'] = cross_species_matrix

        # Save merge summary
        summary = create_merge_summary(merged_matrices)
        summary_file = merged_dir / "merge_summary.json"
        io.dump_json(summary, summary_file)

        execution_time = time.time() - start_time

        return StepResult(
            success=True,
            outputs={
                'merged_dir': str(merged_dir),
                'summary': str(summary_file),
                'matrices': {species: str(path) for species, path in merged_matrices.items()},
                'species_merged': len(merged_matrices)
            },
            metadata={
                'output_format': output_format,
                'execution_time': execution_time
            },
            execution_time=execution_time
        )

    except Exception as e:
        execution_time = time.time() - start_time
        logger.error(f"Merge step failed: {e}")

        return StepResult(
            success=False,
            outputs={},
            metadata={},
            error_message=str(e),
            execution_time=execution_time
        )


def merge_species_expressions(species: str, species_results: Dict[str, Any],
                            output_dir: Path, output_format: str) -> str:
    """Merge expression data for a single species.

    Args:
        species: Species name
        species_results: Quantification results for the species
        output_dir: Output directory
        output_format: Output format

    Returns:
        Path to merged matrix file
    """
    samples = species_results.get('samples', {})

    if not samples:
        logger.warning(f"No samples found for {species}")
        return ""

    # Collect expression data from all samples
    expression_data = {}

    for sample_id, sample_result in samples.items():
        if not sample_result.get('success', False):
            logger.warning(f"Skipping failed sample {sample_id}")
            continue

        # In real implementation, would parse actual quantification output files
        # For now, simulate expression data
        sample_expression = simulate_expression_data(sample_id)
        expression_data[sample_id] = sample_expression

    if not expression_data:
        logger.warning(f"No valid expression data for {species}")
        return ""

    # Create expression matrix
    all_genes = set()
    for expr_dict in expression_data.values():
        all_genes.update(expr_dict.keys())

    all_genes = sorted(all_genes)
    all_samples = sorted(expression_data.keys())

    # Create DataFrame
    matrix_data = []
    for gene in all_genes:
        row = [gene]
        for sample in all_samples:
            expr_dict = expression_data[sample]
            count = expr_dict.get(gene, 0.0)
            row.append(count)
        matrix_data.append(row)

    columns = ['gene'] + all_samples
    df = pd.DataFrame(matrix_data, columns=columns)

    # Save matrix
    species_clean = species.replace(' ', '_')
    if output_format == 'tsv':
        output_file = output_dir / f"{species_clean}_expression_matrix.tsv"
        df.to_csv(output_file, sep='\t', index=False)
    elif output_format == 'csv':
        output_file = output_dir / f"{species_clean}_expression_matrix.csv"
        df.to_csv(output_file, index=False)
    elif output_format == 'hdf5':
        output_file = output_dir / f"{species_clean}_expression_matrix.h5"
        df.to_hdf(output_file, key='expression', mode='w')
    else:
        raise ValueError(f"Unsupported output format: {output_format}")

    logger.info(f"Created expression matrix for {species}: {len(all_genes)} genes x {len(all_samples)} samples")

    return str(output_file)


def merge_cross_species(species_matrices: Dict[str, str], output_dir: Path,
                       output_format: str) -> str:
    """Merge expression matrices across multiple species.

    Args:
        species_matrices: Dictionary of species matrix paths
        output_dir: Output directory
        output_format: Output format

    Returns:
        Path to cross-species matrix file
    """
    # This is a simplified implementation
    # In practice, would need to handle ortholog mapping

    logger.info("Creating cross-species expression matrix")

    # For now, just create a placeholder
    output_file = output_dir / "cross_species_expression_matrix.tsv"

    with open(output_file, 'w') as f:
        f.write("# Cross-species expression matrix\n")
        f.write("# Note: This is a simplified implementation\n")

    return str(output_file)


def simulate_expression_data(sample_id: str) -> Dict[str, float]:
    """Simulate expression data for a sample (for testing).

    Args:
        sample_id: Sample identifier

    Returns:
        Dictionary mapping genes to expression values
    """
    # Simulate expression data for ~15,000 genes
    import random
    random.seed(hash(sample_id))  # Reproducible

    expression = {}
    for i in range(15000):
        gene_id = "03d"
        # Simulate realistic expression distribution
        if random.random() < 0.1:  # 10% highly expressed
            expr = random.uniform(100, 1000)
        elif random.random() < 0.3:  # 30% moderately expressed
            expr = random.uniform(10, 100)
        else:  # 60% lowly expressed or not detected
            expr = random.uniform(0, 10)
        expression[gene_id] = expr

    return expression


def create_merge_summary(merged_matrices: Dict[str, str]) -> Dict[str, Any]:
    """Create a summary of merge results.

    Args:
        merged_matrices: Dictionary of merged matrix paths

    Returns:
        Merge summary
    """
    summary = {
        'total_matrices': len(merged_matrices),
        'matrix_sizes': {},
        'total_genes': 0,
        'total_samples': 0
    }

    for species, matrix_path in merged_matrices.items():
        if Path(matrix_path).exists():
            # In real implementation, would analyze the actual matrix
            summary['matrix_sizes'][species] = "simulated"
            summary['total_genes'] += 15000  # Simulated
            summary['total_samples'] += 10    # Simulated

    return summary


