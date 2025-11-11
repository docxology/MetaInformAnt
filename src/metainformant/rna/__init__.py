"""RNA domain functionality.

This module exposes a thin, modular wrapper around the external `amalgkit`
CLI to support transcriptomic meta-analysis workflows.
"""

from .amalgkit import (
    AmalgkitParams,
    build_amalgkit_command,
    build_cli_args,
    check_cli_available,
    config,
    csca,
    cstmm,
    curate,
    ensure_cli_available,
    getfastq,
    integrate,
    merge,
    metadata,
    quant,
    run_amalgkit,
    sanity,
    select,
)

# Lazy workflow imports to avoid import-time failures on older Python during
# light-weight uses that only need the thin CLI wrappers. Full workflow
# functionality requires Python 3.11+ per project configuration.
try:  # pragma: no cover - exercised in integration tests under py311
    from .workflow import AmalgkitWorkflowConfig, execute_workflow, plan_workflow, load_workflow_config

    _HAS_WORKFLOW = True
except Exception:  # pragma: no cover - defensive for environments <3.11
    _HAS_WORKFLOW = False

# Step functions from steps module
from .steps import (
    # Individual step runners
    run_config,
    run_csca,
    run_cstmm,
    run_curate,
    run_getfastq,
    run_integrate,
    run_merge,
    run_metadata,
    run_quant,
    run_sanity,
    run_select,
    # Processing functions
    run_download_quant_workflow,
    quantify_sample,
    convert_sra_to_fastq,
    delete_sample_fastqs,
)

# Genome preparation functions
from .genome_prep import (
    build_kallisto_index,
    download_cds_fasta_from_ftp,
    download_rna_fasta_from_ftp,
    extract_transcripts_from_gff,
    find_rna_fasta_in_genome_dir,
    get_expected_index_path,
    orchestrate_genome_setup,
    prepare_genome_for_quantification,
    prepare_transcriptome_for_kallisto,
    verify_genome_status,
)

# Monitoring and environment functions (always available)
from .environment import (
    check_amalgkit,
    check_dependencies,
    check_kallisto,
    check_metainformant,
    check_rscript,
    check_sra_toolkit,
    check_virtual_env,
    validate_environment,
)
from .monitoring import (
    analyze_species_status,
    assess_all_species_progress,
    check_active_downloads,
    check_workflow_progress,
    count_quantified_samples,
    find_unquantified_samples,
    get_sample_status,
    initialize_progress_tracking,
)

# Orchestration, cleanup, and discovery functions
from .orchestration import (
    cleanup_unquantified_samples,
    discover_species_configs,
    monitor_workflows,
    run_workflow_for_species,
)
from .cleanup import (
    cleanup_partial_downloads,
    fix_abundance_naming,
    fix_abundance_naming_for_species,
)
from .discovery import (
    generate_config_yaml,
    get_genome_info,
    search_species_with_rnaseq,
)

__all__ = [
    # Amalgkit CLI wrapper
    "AmalgkitParams",
    "build_cli_args",
    "build_amalgkit_command",
    "check_cli_available",
    "ensure_cli_available",
    "run_amalgkit",
    # Amalgkit step functions (high-level)
    "metadata",
    "integrate",
    "config",
    "select",
    "getfastq",
    "quant",
    "merge",
    "cstmm",
    "curate",
    "csca",
    "sanity",
    # Step runners (low-level)
    "run_metadata",
    "run_integrate",
    "run_config",
    "run_select",
    "run_getfastq",
    "run_quant",
    "run_merge",
    "run_cstmm",
    "run_curate",
    "run_csca",
    "run_sanity",
    # Processing functions
    "run_download_quant_workflow",
    "quantify_sample",
    "convert_sra_to_fastq",
    "delete_sample_fastqs",
    # Genome preparation functions
    "build_kallisto_index",
    "download_cds_fasta_from_ftp",
    "download_rna_fasta_from_ftp",
    "extract_transcripts_from_gff",
    "find_rna_fasta_in_genome_dir",
    "get_expected_index_path",
    "orchestrate_genome_setup",
    "prepare_genome_for_quantification",
    "prepare_transcriptome_for_kallisto",
    "verify_genome_status",
    # Monitoring
    "analyze_species_status",
    "assess_all_species_progress",
    "check_active_downloads",
    "check_workflow_progress",
    "count_quantified_samples",
    "find_unquantified_samples",
    "get_sample_status",
    "initialize_progress_tracking",
    # Orchestration
    "cleanup_unquantified_samples",
    "discover_species_configs",
    "monitor_workflows",
    "run_workflow_for_species",
    # Cleanup
    "cleanup_partial_downloads",
    "fix_abundance_naming",
    "fix_abundance_naming_for_species",
    # Discovery
    "generate_config_yaml",
    "get_genome_info",
    "search_species_with_rnaseq",
    # Environment
    "check_amalgkit",
    "check_dependencies",
    "check_kallisto",
    "check_metainformant",
    "check_rscript",
    "check_sra_toolkit",
    "check_virtual_env",
    "validate_environment",
]

if _HAS_WORKFLOW:
    __all__ += [
        "AmalgkitWorkflowConfig",
        "plan_workflow",
        "execute_workflow",
        "load_workflow_config",
    ]
