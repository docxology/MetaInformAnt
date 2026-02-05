"""GWAS data retrieval and download utilities.

This module provides functions for downloading reference genomes, variant databases,
SRA sequencing data, and genome annotation for GWAS analysis.
"""

from metainformant.gwas.data.download import (
    download_annotation,
    download_reference_genome,
    download_variant_data,
    download_variant_database,
)
from metainformant.gwas.data.genome import (
    AMEL_HAV3_CHROM_SIZES,
    AMEL_HAV3_CHROMOSOMES,
    get_chromosome_order,
    normalize_chromosome_name,
    parse_gff3_genes,
)
from metainformant.gwas.data.metadata import (
    get_geographic_coordinates,
    get_population_labels,
    load_sample_metadata,
    merge_metadata_with_phenotypes,
    validate_metadata,
)
from metainformant.gwas.data.sra_download import (
    batch_download_sra,
    check_sra_tools_available,
    download_sra_biosample,
    download_sra_experiment,
    download_sra_project,
    download_sra_run,
    download_sra_with_retry,
    find_sra_data_by_phenotype,
    prefetch_sra_metadata,
    search_sra_for_organism,
    validate_sra_download,
)

__all__ = [
    # Reference data downloads
    "download_reference_genome",
    "download_variant_database",
    "download_annotation",
    "download_variant_data",
    # Genome mapping
    "AMEL_HAV3_CHROMOSOMES",
    "AMEL_HAV3_CHROM_SIZES",
    "normalize_chromosome_name",
    "get_chromosome_order",
    "parse_gff3_genes",
    # SRA downloads
    "download_sra_experiment",
    "download_sra_biosample",
    "download_sra_project",
    "download_sra_run",
    "download_sra_with_retry",
    "batch_download_sra",
    # SRA search and utilities
    "search_sra_for_organism",
    "find_sra_data_by_phenotype",
    "prefetch_sra_metadata",
    "validate_sra_download",
    "check_sra_tools_available",
    # Sample metadata
    "load_sample_metadata",
    "merge_metadata_with_phenotypes",
    "validate_metadata",
    "get_population_labels",
    "get_geographic_coordinates",
]
