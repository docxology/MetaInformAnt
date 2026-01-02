"""GWAS workflow orchestration and configuration management.

This module provides high-level functions for managing GWAS analysis workflows,
including configuration loading, workflow planning, and execution orchestration.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


class GWASWorkflowConfig:
    """Configuration class for GWAS workflows."""

    def __init__(self,
                 work_dir: Path,
                 threads: int = 8,
                 vcf_path: Optional[str] = None,
                 phenotype_path: Optional[str] = None,
                 **kwargs):
        """Initialize GWAS workflow configuration.

        Args:
            work_dir: Working directory for GWAS analysis
            threads: Number of threads to use
            vcf_path: Path to VCF file with genotypes
            phenotype_path: Path to phenotype file
            **kwargs: Additional configuration options
        """
        self.work_dir = Path(work_dir)
        self.threads = threads
        self.vcf_path = vcf_path
        self.phenotype_path = phenotype_path
        self.extra_config = kwargs

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'work_dir': str(self.work_dir),
            'threads': self.threads,
            'vcf_path': self.vcf_path,
            'phenotype_path': self.phenotype_path,
            **self.extra_config
        }

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> GWASWorkflowConfig:
        """Create configuration from dictionary."""
        return cls(**config_dict)


def load_gwas_config(config_file: Union[str, Path]) -> Dict[str, Any]:
    """Load GWAS configuration from YAML/TOML/JSON file.

    Args:
        config_file: Path to configuration file

    Returns:
        Configuration dictionary

    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If config format is invalid
    """
    config_path = Path(config_file)

    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    try:
        if config_path.suffix.lower() in ['.yaml', '.yml']:
            config = io.load_yaml(str(config_path))
        elif config_path.suffix.lower() == '.toml':
            config = io.load_toml(str(config_path))
        elif config_path.suffix.lower() == '.json':
            config = io.load_json(str(config_path))
        else:
            raise ValueError(f"Unsupported config format: {config_path.suffix}")

        return config

    except Exception as e:
        raise ValueError(f"Error loading configuration from {config_path}: {e}")


def execute_gwas_workflow(config: Dict[str, Any], *,
                         check: bool = False) -> Dict[str, Any]:
    """Execute the complete GWAS workflow.

    Args:
        config: GWAS configuration dictionary
        check: If True, only validate configuration without executing

    Returns:
        Workflow results dictionary
    """
    logger.info("Starting GWAS workflow execution")

    if check:
        logger.info("Running in check mode - validating configuration only")
        # Validate configuration
        errors = validate_gwas_config(config)
        if errors:
            return {
                'success': False,
                'errors': errors,
                'config_valid': False
            }
        return {
            'success': True,
            'errors': [],
            'config_valid': True
        }

    # Execute workflow steps
    results = {
        'success': True,
        'steps_completed': [],
        'errors': [],
        'outputs': {}
    }

    try:
        # Step 1: Load and validate data
        logger.info("Step 1: Loading and validating data")
        vcf_data = parse_vcf_full(config['vcf_path'])
        results['steps_completed'].append('data_loading')

        # Step 2: Apply QC filters
        logger.info("Step 2: Applying quality control filters")
        qc_config = config.get('quality_control', {})
        filtered_data = apply_qc_filters(
            vcf_data,
            min_maf=qc_config.get('min_maf', 0.01),
            max_missing=qc_config.get('max_missing', 0.1),
            min_hwe_p=qc_config.get('min_hwe_p', 1e-6)
        )
        results['steps_completed'].append('quality_control')

        # Step 3: Population structure analysis
        logger.info("Step 3: Analyzing population structure")
        genotype_matrix = _extract_genotype_matrix(filtered_data)
        pca_result = compute_pca(genotype_matrix)
        kinship_matrix = compute_kinship_matrix(genotype_matrix)
        results['steps_completed'].append('population_structure')

        # Step 4: Association testing
        logger.info("Step 4: Performing association testing")
        phenotypes = _load_phenotypes(config['phenotype_path'])

        association_results = []
        for i, phenotype in enumerate(phenotypes):
            result = association_test_linear(
                filtered_data['genotypes'][:, i],  # Use first trait for now
                phenotype,
                covariates=pca_result[0][:, :10]  # Use first 10 PCs as covariates
            )
            association_results.append(result)

        results['steps_completed'].append('association_testing')
        results['outputs']['association_results'] = association_results

        # Step 5: Multiple testing correction
        logger.info("Step 5: Applying multiple testing correction")
        p_values = [r['p_value'] for r in association_results]
        corrected_results = fdr_correction(p_values)
        results['steps_completed'].append('multiple_testing_correction')

        # Step 6: Generate plots
        logger.info("Step 6: Generating visualization plots")
        plot_results = generate_all_plots(
            association_results,
            config.get('output_dir', config['work_dir']),
            pca_file=Path(config['work_dir']) / 'pca_results.npy',
            kinship_file=Path(config['work_dir']) / 'kinship_matrix.npy',
            vcf_file=Path(config['vcf_path'])
        )
        results['steps_completed'].append('visualization')
        results['outputs']['plots'] = plot_results

        logger.info("GWAS workflow completed successfully")

    except Exception as e:
        logger.error(f"GWAS workflow failed: {e}")
        results['success'] = False
        results['errors'].append(str(e))

    return results


def validate_gwas_config(config: Dict[str, Any]) -> List[str]:
    """Validate GWAS configuration.

    Args:
        config: Configuration dictionary to validate

    Returns:
        List of validation error messages
    """
    errors = []

    # Required fields
    required_fields = ['work_dir', 'vcf_path', 'phenotype_path']
    for field in required_fields:
        if field not in config:
            errors.append(f"Missing required field: {field}")

    # Check file paths
    if 'vcf_path' in config and not Path(config['vcf_path']).exists():
        errors.append(f"VCF file not found: {config['vcf_path']}")

    if 'phenotype_path' in config and not Path(config['phenotype_path']).exists():
        errors.append(f"Phenotype file not found: {config['phenotype_path']}")

    # Check work directory
    if 'work_dir' in config:
        work_dir = Path(config['work_dir'])
        if work_dir.exists() and not work_dir.is_dir():
            errors.append(f"work_dir exists but is not a directory: {work_dir}")

    return errors


def _extract_genotype_matrix(vcf_data: Dict[str, Any]) -> List[List[int]]:
    """Extract genotype matrix from VCF data.

    Args:
        vcf_data: Parsed VCF data from parse_vcf_full()

    Returns:
        Genotype matrix (variants x samples)

    Raises:
        ValueError: If VCF data is malformed or missing genotypes
    """
    if 'genotypes' not in vcf_data:
        raise ValueError("VCF data missing genotypes field")

    genotypes = vcf_data['genotypes']
    if not genotypes:
        raise ValueError("No genotypes found in VCF data")

    # Ensure all genotype rows have the same length (same number of samples)
    sample_count = len(vcf_data.get('samples', []))
    if sample_count == 0:
        raise ValueError("No samples found in VCF data")

    for i, variant_genotypes in enumerate(genotypes):
        if len(variant_genotypes) != sample_count:
            raise ValueError(f"Variant {i} has {len(variant_genotypes)} genotypes but {sample_count} samples expected")

    logger.info(f"Extracted genotype matrix: {len(genotypes)} variants x {sample_count} samples")
    return genotypes


def _load_phenotypes(phenotype_path: str | Path) -> List[float]:
    """Load phenotype data from file.

    Supports common phenotype file formats:
    - CSV/TSV files with headers (looks for 'phenotype', 'trait', or column 1)
    - Plain text files with one value per line
    - Files with sample_id,phenotype format

    Args:
        phenotype_path: Path to phenotype file

    Returns:
        List of phenotype values (one per sample)

    Raises:
        FileNotFoundError: If phenotype file doesn't exist
        ValueError: If phenotype file format is invalid or empty
    """
    phenotype_path = Path(phenotype_path)
    if not phenotype_path.exists():
        raise FileNotFoundError(f"Phenotype file not found: {phenotype_path}")

    phenotypes = []
    logger.info(f"Loading phenotypes from: {phenotype_path}")

    try:
        # Read the file
        with open(phenotype_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        if not lines:
            raise ValueError("Phenotype file is empty")

        # Try to detect format
        first_line = lines[0]

        # Check if it's a CSV/TSV format
        if ',' in first_line or '\t' in first_line:
            delimiter = ',' if ',' in first_line else '\t'
            parts = first_line.split(delimiter)

            # Check if first line looks like header
            # Only consider it a header if it has multiple columns and first column looks like an ID
            has_header = False
            if len(parts) >= 2:
                first_col = parts[0].strip()
                second_col = parts[1].strip()

                # Check if first column looks like a sample ID (contains letters) and second is numeric
                if (any(c.isalpha() for c in first_col) and
                    second_col.replace('.', '').replace('-', '').isdigit()):
                    # This looks like sample_id,phenotype format, not a header
                    has_header = False
                elif any(keyword in first_col.lower() for keyword in ['sample', 'id', 'phenotype', 'trait', 'subject', 'individual']):
                    # First column contains header-like keywords
                    has_header = True

            pheno_col_idx = -1
            if has_header:
                # Look for phenotype column in header
                header = [col.strip().lower() for col in parts]
                for i, col in enumerate(header):
                    if col in ['phenotype', 'trait', 'pheno', 'phenotypes']:
                        pheno_col_idx = i
                        break

                # If no phenotype column found, assume second column (after sample ID)
                if pheno_col_idx == -1 and len(parts) >= 2:
                    pheno_col_idx = 1

                # Parse data rows (skip header)
                data_lines = lines[1:]
            else:
                # No header, assume phenotype is in second column
                pheno_col_idx = 1
                data_lines = lines

            if pheno_col_idx == -1:
                raise ValueError("Could not identify phenotype column")

            # Parse data rows
            for line in data_lines:
                if line.strip():
                    parts = line.split(delimiter)
                    if len(parts) > pheno_col_idx:
                        try:
                            value = float(parts[pheno_col_idx].strip())
                            phenotypes.append(value)
                        except ValueError:
                            logger.warning(f"Could not parse phenotype value: {parts[pheno_col_idx]}")
                            continue

        else:
            # Assume simple format: one phenotype per line
            for line in lines:
                try:
                    value = float(line.strip())
                    phenotypes.append(value)
                except ValueError:
                    logger.warning(f"Could not parse phenotype value: {line}")
                    continue

        if not phenotypes:
            raise ValueError("No valid phenotype values found in file")

        logger.info(f"Loaded {len(phenotypes)} phenotypes")
        return phenotypes

    except Exception as e:
        raise ValueError(f"Error loading phenotypes from {phenotype_path}: {e}")


# Import required functions from GWAS modules
from .quality import parse_vcf_full, apply_qc_filters
from .structure import compute_pca, compute_kinship_matrix
from .association import association_test_linear
from .correction import fdr_correction
from .visualization import generate_all_plots


def run_gwas(vcf_path: Union[str, Path], phenotype_path: Union[str, Path],
             config: Dict[str, Any], output_dir: Union[str, Path] | None = None) -> Dict[str, Any]:
    """Run complete GWAS workflow.

    Args:
        vcf_path: Path to VCF file
        phenotype_path: Path to phenotype file
        config: GWAS configuration dictionary
        output_dir: Output directory (optional)

    Returns:
        Dictionary with GWAS results and metadata

    Raises:
        ValueError: If required parameters are missing
        FileNotFoundError: If input files don't exist
    """
    import tempfile
    from pathlib import Path

    logger.info("Starting GWAS workflow")

    # Validate inputs
    vcf_path = Path(vcf_path)
    phenotype_path = Path(phenotype_path)

    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    if not phenotype_path.exists():
        raise FileNotFoundError(f"Phenotype file not found: {phenotype_path}")

    # Create output directory
    if output_dir is None:
        temp_dir = tempfile.mkdtemp()
        output_dir = Path(temp_dir)
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    results = {
        'config': config,
        'output_dir': str(output_dir),
        'steps_completed': [],
        'results': {}
    }

    try:
        # Step 1: Load and parse VCF
        logger.info("Parsing VCF file")
        vcf_data = parse_vcf_full(vcf_path)
        results['steps_completed'].append('parse_vcf')
        results['results']['vcf_summary'] = {
            'num_variants': len(vcf_data.get('variants', [])),
            'num_samples': len(vcf_data.get('samples', []))
        }

        # Step 2: Apply QC filters
        logger.info("Applying QC filters")
        qc_config = config.get('quality_control', {})
        filtered_data = apply_qc_filters(
            vcf_data,
            min_maf=qc_config.get('min_maf', 0.01),
            max_missing=qc_config.get('max_missing', 0.1),
            min_hwe_p=qc_config.get('min_hwe_p', 1e-6)
        )
        results['steps_completed'].append('qc_filters')
        results['results']['qc_summary'] = {
            'variants_after_qc': len(filtered_data.get('variants', [])),
            'samples_after_qc': len(filtered_data.get('samples', []))
        }

        # Step 3: Population structure analysis
        logger.info("Computing population structure")
        genotypes = filtered_data.get('genotypes', [])
        if genotypes:
            pca_result = compute_pca(genotypes, n_components=min(10, len(genotypes)))
            kinship_matrix = compute_kinship_matrix(genotypes)
            results['steps_completed'].append('population_structure')
            results['results']['pca'] = pca_result
            results['results']['kinship'] = kinship_matrix

        # Step 4: Load phenotypes
        logger.info("Loading phenotypes")
        phenotypes = load_phenotypes(phenotype_path)
        results['steps_completed'].append('load_phenotypes')
        results['results']['phenotype_summary'] = {
            'num_samples': len(phenotypes),
            'trait_name': config.get('trait_name', 'unknown')
        }

        # Step 5: Association testing
        logger.info("Running association tests")
        assoc_results = []
        if genotypes and phenotypes:
            for i, genotype in enumerate(genotypes):
                try:
                    result = association_test_linear(genotype, phenotypes)
                    result['variant_id'] = f'variant_{i}'
                    assoc_results.append(result)
                except Exception as e:
                    logger.warning(f"Association test failed for variant {i}: {e}")

        results['steps_completed'].append('association_testing')
        results['results']['association_results'] = assoc_results

        # Step 6: Multiple testing correction
        logger.info("Applying multiple testing correction")
        if assoc_results:
            p_values = [r.get('p_value', 1.0) for r in assoc_results]
            _, q_values = fdr_correction(p_values)

            for i, q_val in enumerate(q_values):
                if i < len(assoc_results):
                    assoc_results[i]['q_value'] = q_val

        results['steps_completed'].append('multiple_testing_correction')

        # Step 7: Generate plots
        logger.info("Generating visualization plots")
        try:
            plot_results = generate_all_plots(
                association_results=str(output_dir / "association_results.json"),
                output_dir=output_dir,
                significance_threshold=config.get('significance_threshold', 5e-8)
            )
            results['steps_completed'].append('visualization')
            results['results']['plots'] = plot_results
        except Exception as e:
            logger.warning(f"Plot generation failed: {e}")

        results['status'] = 'completed'
        logger.info("GWAS workflow completed successfully")

    except Exception as e:
        results['status'] = 'failed'
        results['error'] = str(e)
        logger.error(f"GWAS workflow failed: {e}")
        raise

    return results
