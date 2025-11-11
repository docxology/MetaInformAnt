from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping

from ..core.config import apply_env_overrides, load_mapping_from_file
from ..core.io import dump_json, ensure_directory, read_delimited, write_delimited, write_jsonl
from ..core.errors import error_context
from ..core.logging import get_logger, log_with_metadata
from ..core.paths import expand_and_resolve
from . import steps as _steps_mod
from .amalgkit import AmalgkitParams, build_amalgkit_command, check_cli_available, ensure_cli_available
from .deps import check_step_dependencies

logger = get_logger(__name__)


@dataclass
class AmalgkitWorkflowConfig:
    """Configuration for amalgkit RNA-seq workflow execution.
    
    Attributes:
        work_dir: Base working directory for all workflow outputs
        threads: Number of threads to use (default: 6)
        species_list: List of species names to process
        log_dir: Optional directory for log files
        manifest_path: Optional path to workflow manifest file
        per_step: Per-step parameter overrides
        auto_install_amalgkit: Whether to attempt automatic installation
        genome: Genome configuration dictionary
        filters: Filter criteria for sample selection
    """
    work_dir: Path
    threads: int = 6
    species_list: list[str] = field(default_factory=list)
    # Additional common parameters can be added as needed
    log_dir: Path | None = None
    manifest_path: Path | None = None
    per_step: dict[str, AmalgkitParams] = field(default_factory=dict)
    auto_install_amalgkit: bool = False
    genome: dict[str, Any] | None = (
        None  # e.g., { accession: GCF_*, index_out: path, dest_dir: path, ftp_url: str, include: [...] }
    )
    filters: dict[str, Any] = field(default_factory=dict)  # e.g., { require_tissue: true }

    def to_common_params(self) -> AmalgkitParams:
        """Convert configuration to common parameters for all steps.
        
        Returns:
            Dictionary of common parameters (threads, species-list)
        """
        params: dict[str, Any] = {"threads": self.threads}
        if self.species_list:
            params["species-list"] = list(self.species_list)
        return params


def _find_repo_root(config_file: Path) -> Path:
    """Find repository root directory by walking up from config file location.
    
    Looks for common repository markers: .git, pyproject.toml, or setup.py.
    If not found, uses the directory containing the config file as fallback.
    
    Args:
        config_file: Path to config file
        
    Returns:
        Path to repository root directory
    """
    config_path = Path(config_file).resolve()
    current = config_path.parent
    
    # Walk up directory tree looking for repo markers
    markers = [".git", "pyproject.toml", "setup.py", ".cursorrules"]
    for _ in range(10):  # Limit search to 10 levels up
        for marker in markers:
            if (current / marker).exists():
                return current
        if current.parent == current:  # Reached filesystem root
            break
        current = current.parent
    
    # Fallback: use config directory's parent (assuming config/ is in repo root)
    if "config" in config_path.parts:
        config_index = config_path.parts.index("config")
        return Path(*config_path.parts[:config_index])
    
    # Last resort: use current working directory
    return Path.cwd()


def _resolve_path_relative_to_repo(path_str: str | Path, repo_root: Path) -> Path:
    """Resolve a path relative to repository root.
    
    If path is absolute, returns it as-is (after expand_and_resolve).
    If path is relative, resolves it relative to repo_root.
    
    Args:
        path_str: Path string (can be relative or absolute)
        repo_root: Repository root directory
        
    Returns:
        Resolved absolute Path
    """
    path = Path(path_str)
    
    # If absolute, just expand and resolve
    if path.is_absolute():
        return expand_and_resolve(path)
    
    # If relative, resolve relative to repo root
    return expand_and_resolve(repo_root / path)


def _build_default_search_string(species_list: list[str]) -> str | None:
    if not species_list:
        return None
    if len(species_list) == 1:
        sp = species_list[0].replace("_", " ")
        return f'"{sp}"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]'
    # Multiple species: join with OR
    parts = [f'"{sp.replace("_", " ")}"[Organism]' for sp in species_list]
    return f'({" OR ".join(parts)}) AND RNA-Seq[Strategy] AND Illumina[Platform]'


def _apply_step_defaults(config: AmalgkitWorkflowConfig) -> None:
    """Fill in sensible defaults for amalgkit steps when not provided in config."""
    ps = config.per_step
    # metadata defaults
    md = dict(ps.get("metadata", {}))
    md.setdefault("out_dir", str(config.work_dir))
    search = _build_default_search_string(config.species_list)
    if search and not md.get("search_string"):
        md["search_string"] = search
    email = os.environ.get("NCBI_EMAIL")
    if email and not md.get("entrez_email"):
        md["entrez_email"] = email
    ps["metadata"] = md

    # directories
    defaults = {
        "integrate": {"out_dir": str(config.work_dir), "fastq_dir": str(config.work_dir / "fastq")},
        "config": {"out_dir": str(config.work_dir), "config": "base"},
        "select": {"out_dir": str(config.work_dir), "config_dir": str(config.work_dir / "config_base")},
        "getfastq": {"out_dir": str(config.work_dir / "fastq")},
        "quant": {"out_dir": str(config.work_dir / "quant")},
        "merge": {"out_dir": str(config.work_dir.parent / "merged")},
        "cstmm": {"out_dir": str(config.work_dir / "cstmm")},
        "curate": {"out_dir": str(config.work_dir / "curate")},
        "csca": {"out_dir": str(config.work_dir / "csca")},
        "sanity": {"out_dir": str(config.work_dir)},
    }
    for step, d in defaults.items():
        cur = dict(ps.get(step, {}))
        for k, v in d.items():
            cur.setdefault(k, v)
        ps[step] = cur


def plan_workflow(config: AmalgkitWorkflowConfig) -> list[tuple[str, AmalgkitParams]]:
    """Return an ordered list of (subcommand, params) representing a full run.

    This does not execute anything; it allows dry inspection and TDD.
    
    Note: integrate runs AFTER getfastq to integrate downloaded FASTQs into metadata.
    """
    common = config.to_common_params()

    def merge_params(extra: Mapping[str, Any] | None = None) -> AmalgkitParams:
        """Merge common parameters with step-specific overrides.
        
        Args:
            extra: Step-specific parameter overrides
            
        Returns:
            Merged parameter dictionary
        """
        if not extra:
            return dict(common)
        merged = dict(common)
        merged.update(extra)
        return merged

    steps: list[tuple[str, AmalgkitParams]] = []
    ordered = [
        "metadata",
        "config",
        "select",
        "getfastq",
        "integrate",  # Moved after getfastq to integrate downloaded FASTQs
        "quant",
        "merge",
        "cstmm",
        "curate",
        "csca",
        "sanity",
    ]
    for name in ordered:
        extra = config.per_step.get(name, {})
        steps.append((name, merge_params(extra)))
    return steps


def plan_workflow_with_params(
    config: AmalgkitWorkflowConfig,
    step_params: dict[str, AmalgkitParams],
) -> list[tuple[str, AmalgkitParams]]:
    """Plan workflow while merging per-step params on top of common ones."""
    common = config.to_common_params()

    ordered = [
        "metadata",
        "config",
        "select",
        "getfastq",
        "integrate",  # Moved after getfastq
        "quant",
        "merge",
        "cstmm",
        "curate",
        "csca",
        "sanity",
    ]
    steps: list[tuple[str, AmalgkitParams]] = []
    for name in ordered:
        params = dict(common)
        params.update(step_params.get(name, {}))
        steps.append((name, params))
    return steps


def _write_manifest_records(path: Path, records: list[dict[str, Any]]) -> None:
    ensure_directory(path.parent)
    write_jsonl(records, path)


def _sanitize_params_for_subcommand(subcommand: str, params: Mapping[str, Any]) -> dict[str, Any]:
    """Drop unsupported flags for specific subcommands.

    Currently, `amalgkit metadata` accepts: out_dir, redo, search_string, entrez_email, resolve_names.
    `amalgkit select` accepts: out_dir, metadata, sample_group, config_dir, mark_missing_rank, min_nspots, max_sample, mark_redundant_biosamples.
    """
    allowed_map: dict[str, set[str]] = {
        "metadata": {"out_dir", "redo", "search_string", "entrez_email", "resolve_names"},
        "integrate": {"out_dir", "metadata", "threads", "fastq_dir", "id", "id_list"},
        "config": {"out_dir", "config"},
        "select": {"out_dir", "metadata", "sample_group", "config_dir", "mark_missing_rank", "min_nspots", "max_sample", "mark_redundant_biosamples"},
        # getfastq: allow robust flags pass-through
        "getfastq": {
            "out_dir",
            "metadata",
            "threads",
            "redo",
            "entrez_email",
            "id",
            "id_list",
            # robustness/diagnostics
            "pfd",
            "pfd_exe",
            "prefetch_exe",
            "fastp",
            "fastp_exe",
            "pfd_print",
            "fastp_print",
            # cloud source toggles
            "aws",
            "ncbi",
            "gcp",
            # accelerate is internal to our step runner, not passed to amalgkit CLI
        },
        "quant": {"out_dir", "metadata", "threads", "redo", "index_dir"},
        "merge": {"out", "out_dir", "metadata"},
        "cstmm": {"out_dir", "metadata", "orthogroup_table"},
        "curate": {"out_dir", "metadata", "input_dir", "sample_group", "sample_group_color", "redo"},
        "csca": {"out_dir", "metadata", "sample_group", "sample_group_color", "orthogroup_table"},
        "sanity": {"out_dir", "metadata", "index"},
    }
    allowed = allowed_map.get(subcommand)
    if allowed is None:
        return dict(params)
    return {k: v for k, v in params.items() if k in allowed}


def _validate_genome_ready(config: AmalgkitWorkflowConfig) -> dict[str, Any]:
    """Validate that genome and kallisto index are ready for quantification.
    
    Checks:
    - Genome directory exists and contains transcriptome FASTA
    - Kallisto index exists if build_index is enabled in quant config
    
    Args:
        config: Workflow configuration
        
    Returns:
        Dictionary with validation results:
        - ready: bool - True if genome/index ready
        - genome_exists: bool - True if genome directory exists
        - transcriptome_exists: bool - True if transcriptome FASTA found
        - index_exists: bool - True if kallisto index exists (or not required)
        - genome_dir: Path | None - Path to genome directory
        - transcriptome_path: Path | None - Path to transcriptome FASTA
        - index_path: Path | None - Path to kallisto index
        - errors: list[str] - List of error messages
    """
    result: dict[str, Any] = {
        "ready": False,
        "genome_exists": False,
        "transcriptome_exists": False,
        "index_exists": False,
        "genome_dir": None,
        "transcriptome_path": None,
        "index_path": None,
        "errors": [],
    }
    
    if not config.genome:
        result["ready"] = True  # No genome config means no validation needed
        return result
    
    from .genome_prep import find_rna_fasta_in_genome_dir, get_expected_index_path
    
    # Check genome directory
    acc = str(config.genome.get("accession", ""))
    default_dest = config.work_dir.parent / "genome"
    dest_dir_str = config.genome.get("dest_dir", default_dest)
    dest_dir = Path(dest_dir_str) if isinstance(dest_dir_str, str) else dest_dir_str
    
    if not dest_dir.exists():
        result["errors"].append(f"Genome directory does not exist: {dest_dir}")
        return result
    
    result["genome_exists"] = True
    result["genome_dir"] = dest_dir
    
    # Check for transcriptome FASTA
    transcriptome_path = find_rna_fasta_in_genome_dir(dest_dir, acc)
    if not transcriptome_path:
        result["errors"].append(f"Transcriptome FASTA not found in genome directory: {dest_dir}")
        return result
    
    result["transcriptome_exists"] = True
    result["transcriptome_path"] = transcriptome_path
    
    # Check for kallisto index if build_index is enabled
    quant_params = dict(config.per_step.get("quant", {}))
    build_index = quant_params.get("build_index", False)
    if isinstance(build_index, str):
        build_index = build_index.lower() in {"yes", "true", "1"}
    
    if build_index:
        species_name = config.species_list[0] if config.species_list else "unknown"
        index_path = get_expected_index_path(config.work_dir, species_name)
        
        if index_path.exists():
            result["index_exists"] = True
            result["index_path"] = index_path
            result["ready"] = True
        else:
            result["errors"].append(f"Kallisto index not found (expected: {index_path})")
    else:
        # Index not required
        result["index_exists"] = True
        result["ready"] = True
    
    return result


def execute_workflow(config: AmalgkitWorkflowConfig, *, check: bool = False) -> list[int]:
    """Execute the full amalgkit workflow in order.

    This workflow provides end-to-end functionality:
    1. Automatic genome download → transcriptome preparation → kallisto indexing (if genome config exists)
    2. Metadata retrieval for available samples
    3. Immediate per-sample processing: download → quantify → delete FASTQ (default behavior)
    
    Genome setup is automatic: if genome config exists but genome/index is missing,
    the workflow will automatically download and prepare everything before proceeding.

    Returns a list of return codes per step in order.
    """
    # Log workflow start with metadata
    log_with_metadata(
        logger,
        "Starting RNA workflow execution",
        {
            "work_dir": str(config.work_dir),
            "threads": config.threads,
            "num_species": len(config.species_list) if config.species_list else 0,
            "has_genome_config": config.genome is not None,
        },
    )
    ensure_directory(config.work_dir)
    manifest_path = config.manifest_path or (config.work_dir / "amalgkit.manifest.jsonl")
    logger.info(f"execute_workflow: Manifest path: {manifest_path}")

    manifest_records: list[dict[str, Any]] = []
    return_codes: list[int] = []

    # Ensure defaults are applied for required amalgkit params
    _apply_step_defaults(config)

    # 1) Genome validation and automatic setup FIRST (before any sample processing)
    if config.genome:
        from ..dna.ncbi import download_genome_package_best_effort
        from .genome_prep import prepare_genome_for_quantification
        
        # Validate genome/index status
        validation = _validate_genome_ready(config)
        
        if validation["ready"]:
            logger.info("Genome and index already ready, skipping setup")
            # Use existing paths for quant params
            quant_params = dict(config.per_step.get("quant", {}))
            if validation.get("transcriptome_path"):
                quant_params["fasta_dir"] = str(Path(validation["transcriptome_path"]).parent)
            if validation.get("index_path"):
                quant_params["index_dir"] = str(Path(validation["index_path"]).parent)
            config.per_step["quant"] = quant_params
        else:
            logger.info("Genome/index not ready, running automatic setup...")
            logger.info(f"Validation status: genome_exists={validation['genome_exists']}, "
                       f"transcriptome_exists={validation['transcriptome_exists']}, "
                       f"index_exists={validation['index_exists']}")
            if validation["errors"]:
                logger.info(f"Validation errors: {validation['errors']}")
            
            acc = str(config.genome.get("accession", ""))
            include = config.genome.get("include") or ["gff3", "rna", "cds", "protein", "genome", "seq-report"]
            ftp_url = config.genome.get("ftp_url")
            default_dest = config.work_dir.parent / "genome"
            dest_dir_str = config.genome.get("dest_dir", default_dest)
            # dest_dir is already resolved in load_workflow_config, but handle both string and Path
            dest_dir = Path(dest_dir_str) if isinstance(dest_dir_str, str) else dest_dir_str
            dest_dir = dest_dir.resolve()  # Ensure absolute path
            
            # Step 1: Download genome if missing
            genome_dir_path = None
            if not validation["genome_exists"] or not validation["transcriptome_exists"]:
                logger.info(f"Downloading genome package to {dest_dir}...")
                start_ts = datetime.utcnow()
                dl_rec = download_genome_package_best_effort(acc, dest_dir, include=include, ftp_url=ftp_url)
                end_ts = datetime.utcnow()
                genome_dir_path = Path(dl_rec.get("extracted_dir") or str(dest_dir))
                
                genome_rec = {
                    "step": "genome-download",
                    "return_code": int(dl_rec.get("return_code", 0)),
                    "stdout_bytes": len(dl_rec.get("stdout", "")),
                    "stderr_bytes": len(dl_rec.get("stderr", "")),
                    "started_utc": start_ts.isoformat() + "Z",
                    "finished_utc": end_ts.isoformat() + "Z",
                    "duration_seconds": max(0.0, (end_ts - start_ts).total_seconds()),
                    "work_dir": str(config.work_dir),
                    "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                    "params": {"accession": acc, "include": include, "dest_dir": str(dest_dir), "ftp_url": ftp_url or ""},
                    "command": dl_rec.get("command", dl_rec.get("url", dl_rec.get("method", ""))),
                    "extracted_dir": str(genome_dir_path),
                    "zip_path": dl_rec.get("zip_path", ""),
                    "method": dl_rec.get("method", ""),
                }
                manifest_records.append(genome_rec)
                return_codes.append(genome_rec["return_code"])
            else:
                # Genome already exists, use existing directory
                genome_dir_path = validation["genome_dir"]
                logger.info(f"Genome already exists at {genome_dir_path}, skipping download")
            
            # Step 2: Prepare transcriptome and build index if needed
            if genome_dir_path and (not validation["transcriptome_exists"] or not validation["index_exists"]):
                species_name = config.species_list[0] if config.species_list else "unknown"
                
                # Check if build_index is enabled in quant config
                quant_params = dict(config.per_step.get("quant", {}))
                build_index = quant_params.get("build_index", False)
                if isinstance(build_index, str):
                    build_index = build_index.lower() in {"yes", "true", "1"}
                
                # Only build index if it doesn't exist and is requested
                should_build_index = build_index and not validation["index_exists"]
                
                logger.info(f"Preparing transcriptome (build_index={should_build_index})...")
                prep_start_ts = datetime.utcnow()
                prep_result = prepare_genome_for_quantification(
                    genome_dir_path,
                    species_name,
                    config.work_dir,
                    accession=acc,
                    build_index=should_build_index,
                    kmer_size=31,
                )
                prep_end_ts = datetime.utcnow()
                
                prep_rec = {
                    "step": "transcriptome-prepare",
                    "return_code": 0 if prep_result["success"] else 1,
                    "stdout_bytes": 0,
                    "stderr_bytes": 0,
                    "started_utc": prep_start_ts.isoformat() + "Z",
                    "finished_utc": prep_end_ts.isoformat() + "Z",
                    "duration_seconds": max(0.0, (prep_end_ts - prep_start_ts).total_seconds()),
                    "work_dir": str(config.work_dir),
                    "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                    "params": {
                        "species_name": species_name,
                        "genome_dir": str(genome_dir_path),
                        "build_index": should_build_index,
                    },
                    "fasta_path": prep_result.get("fasta_path"),
                    "index_path": prep_result.get("index_path"),
                    "error": prep_result.get("error"),
                }
                manifest_records.append(prep_rec)
                return_codes.append(prep_rec["return_code"])
                
                # Update quant params with fasta_dir and index_dir if available
                if prep_result.get("fasta_path"):
                    quant_params["fasta_dir"] = str(Path(prep_result["fasta_path"]).parent)
                    logger.info(f"Set fasta_dir for quant: {quant_params['fasta_dir']}")
                
                if prep_result.get("index_path"):
                    quant_params["index_dir"] = str(Path(prep_result["index_path"]).parent)
                    logger.info(f"Set index_dir for quant: {quant_params['index_dir']}")
                
                config.per_step["quant"] = quant_params
            elif validation["ready"]:
                # Everything already exists, just set quant params
                quant_params = dict(config.per_step.get("quant", {}))
                if validation.get("transcriptome_path"):
                    quant_params["fasta_dir"] = str(Path(validation["transcriptome_path"]).parent)
                if validation.get("index_path"):
                    quant_params["index_dir"] = str(Path(validation["index_path"]).parent)
                config.per_step["quant"] = quant_params

    # 2) Ensure amalgkit is available; optionally auto-install
    install_record: dict[str, Any] | None = None
    ok, help_or_msg = check_cli_available()
    if not ok and config.auto_install_amalgkit:
        ok, help_or_msg, install_record = ensure_cli_available(auto_install=True)
        if install_record is not None:
            manifest_records.append(
                {
                    "step": "amalgkit-install",
                    "return_code": install_record.get("return_code", -1),
                    "stdout_bytes": len(install_record.get("stdout", "")),
                    "stderr_bytes": len(install_record.get("stderr", "")),
                    "started_utc": datetime.utcnow().isoformat() + "Z",
                    "finished_utc": datetime.utcnow().isoformat() + "Z",
                    "duration_seconds": 0.0,
                    "work_dir": str(config.work_dir),
                    "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                    "params": {},
                    "command": install_record.get("command", ""),
                    "note": "attempted auto-install of amalgkit",
                }
            )
            return_codes.append(install_record.get("return_code", -1))

    if not ok:
        manifest_records.append(
            {
                "step": "preflight",
                "return_code": 127,
                "stdout_bytes": 0,
                "stderr_bytes": len(help_or_msg or ""),
                "started_utc": datetime.utcnow().isoformat() + "Z",
                "finished_utc": datetime.utcnow().isoformat() + "Z",
                "duration_seconds": 0.0,
                "work_dir": str(config.work_dir),
                "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                "params": {},
                "command": "amalgkit -h",
                "note": "amalgkit not available on PATH",
            }
        )
        return_codes.append(127)
        _write_manifest_records(manifest_path, manifest_records)
        return return_codes

    # 3) Inspect metadata table (if present) to inform downstream decisions
    metadata_table_path = config.work_dir / "metadata" / "metadata.tsv"
    has_metadata_rows = False
    try:
        if metadata_table_path.exists():
            rows = list(read_delimited(metadata_table_path, delimiter="\t"))
            has_metadata_rows = len(rows) > 0
    except Exception:
        has_metadata_rows = False

    # 4) Run amalgkit steps
    logger.info("execute_workflow: Planning workflow steps")
    steps = plan_workflow(config)
    logger.info(f"execute_workflow: Planned {len(steps)} steps: {[s[0] for s in steps]}")
    
    # 4.1) Detect if sequential download-quant-delete processing should be used
    # This prevents disk exhaustion by processing one sample at a time
    step_names = [s[0] for s in steps]
    use_sequential = "getfastq" in step_names and "quant" in step_names
    logger.info(f"execute_workflow: use_sequential={use_sequential}, getfastq={'getfastq' in step_names}, quant={'quant' in step_names}")
    
    if use_sequential:
        logger.info("Detected getfastq + quant in workflow: Using batched processing to manage disk space")
        logger.info("Each sample will be: downloaded → quantified → FASTQ deleted before next sample")
        
        # Find getfastq and quant params
        getfastq_params = None
        quant_params = None
        getfastq_idx = None
        quant_idx = None
        
        for idx, (subcommand, params) in enumerate(steps):
            if subcommand == "getfastq":
                getfastq_params = params
                getfastq_idx = idx
            elif subcommand == "quant":
                quant_params = params
                quant_idx = idx
        
        # Run all steps before getfastq normally
        pre_steps = steps[:getfastq_idx] if getfastq_idx is not None else []
        # Steps after quant (merge, curate, etc.)
        post_steps = steps[quant_idx + 1:] if quant_idx is not None else []
        
        logger.info(f"execute_workflow: Pre-steps: {[s[0] for s in pre_steps]}, Post-steps: {[s[0] for s in post_steps]}")
        
        # Run pre-steps (metadata, config, select, etc.)
        logger.info(f"execute_workflow: Executing {len(pre_steps)} pre-steps...")
        for idx, (subcommand, params) in enumerate(pre_steps):
            logger.info(f"execute_workflow: Pre-step {idx+1}/{len(pre_steps)}: {subcommand}")
            runner = _steps_mod.STEP_RUNNERS.get(subcommand)
            if runner is None:  # pragma: no cover - defensive
                raise KeyError(f"No runner registered for step: {subcommand}")
            
            filtered = _sanitize_params_for_subcommand(subcommand, params)
            
            # Apply the same logic as the main loop for these steps
            if subcommand == "select" and bool(config.filters.get("require_tissue", False)):
                meta_path = filtered.get("metadata")
                if not meta_path or meta_path == "inferred":
                    meta_path = str(config.work_dir / "metadata" / "metadata.tsv")
                meta_path_p = Path(str(meta_path))
                if meta_path_p.exists():
                    out_filtered = meta_path_p.with_name("metadata.filtered.tissue.tsv")
                    try:
                        rows = [
                            row
                            for row in read_delimited(meta_path_p, delimiter="\t")
                            if (row.get("tissue", "").strip() != "") and (row.get("run", "").strip() != "")
                        ]
                        write_delimited(rows, out_filtered, delimiter="\t")
                        filtered["metadata"] = str(out_filtered)
                    except Exception:
                        pass
            
            # Execute pre-step
            start_ts = datetime.utcnow()
            result = runner(
                filtered,
                work_dir=config.work_dir,
                log_dir=(config.log_dir or (config.work_dir / "logs")),
                check=check,
            )
            end_ts = datetime.utcnow()
            duration_s = max(0.0, (end_ts - start_ts).total_seconds())
            return_codes.append(result.returncode)
            manifest_records.append(
                {
                    "step": subcommand,
                    "return_code": result.returncode,
                    "stdout_bytes": len(result.stdout or ""),
                    "stderr_bytes": len(result.stderr or ""),
                    "started_utc": start_ts.isoformat() + "Z",
                    "finished_utc": end_ts.isoformat() + "Z",
                    "duration_seconds": duration_s,
                    "work_dir": str(config.work_dir),
                    "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                    "params": dict(filtered),
                    "command": " ".join(build_amalgkit_command(subcommand, filtered)),
                }
            )
            logger.info(f"execute_workflow: Pre-step {subcommand} completed with code {result.returncode}")
        
        logger.info(f"execute_workflow: All {len(pre_steps)} pre-steps completed")
        
        # Now run batched download + quant
        logger.info(f"execute_workflow: Checking for batched processing: getfastq_params={'exists' if getfastq_params else 'None'}, quant_params={'exists' if quant_params else 'None'}")
        if getfastq_params and quant_params:
            logger.info("execute_workflow: Starting batched processing setup")
            # Determine metadata file to use - must have 'run' column (not a pivot table)
            metadata_paths = [
                config.work_dir / "metadata" / "metadata.filtered.clean.tsv",
                config.work_dir / "metadata" / "metadata.filtered.tissue.tsv",
                config.work_dir / "metadata" / "metadata.tsv",
            ]
            metadata_file = None
            for mp in metadata_paths:
                if mp.exists():
                    # Check if this file has a 'run' column (not a pivot table)
                    try:
                        # Explicitly specify TSV delimiter
                        rows = list(read_delimited(mp, delimiter="\t"))
                        if rows and "run" in rows[0]:
                            metadata_file = mp
                            sample_count = len(rows)
                            logger.info(f"✓ Using metadata file for parallel download workflow: {mp} ({sample_count} samples)")
                            break
                        else:
                            cols = list(rows[0].keys()) if rows else []
                            logger.warning(f"✗ Skipping {mp.name}: no 'run' column. Columns: {cols[:5]}...")
                    except Exception as e:
                        import traceback
                        logger.warning(f"✗ Error checking {mp.name}: {e}")
                        logger.debug(traceback.format_exc())
                        continue
            
            if metadata_file is None:
                logger.error("=" * 80)
                logger.error("No metadata file with 'run' column found for parallel download workflow")
                logger.error(f"Checked {len(metadata_paths)} paths:")
                for p in metadata_paths:
                    logger.error(f"  - {p} (exists: {p.exists()})")
                logger.error("=" * 80)
                return_codes.append(1)
            else:
                # Check if parallel processing is requested
                num_workers = getfastq_params.get("num_download_workers") or getfastq_params.get("parallel_workers") or 8
                
                # Check for max_samples limit (from metadata step or config)
                max_samples = None
                if "metadata" in config.per_step:
                    max_samples = config.per_step["metadata"].get("max_samples")
                if max_samples and isinstance(max_samples, str):
                    try:
                        max_samples = int(max_samples)
                    except ValueError:
                        max_samples = None
                
                # Use unified processing function (handles both sequential and parallel modes)
                if num_workers > 1:
                    logger.info(f"🚀 Using parallel download processing with {num_workers} workers")
                    logger.info(f"   {num_workers} samples download in parallel → quantify sequentially → delete FASTQs")
                    logger.info(f"   This increases throughput while managing disk space")
                else:
                    logger.info(f"🚀 Using immediate per-sample processing")
                    logger.info(f"   Each sample: download → immediately quantify → immediately delete FASTQs → next sample")
                    logger.info(f"   This ensures maximum disk efficiency: only one sample's FASTQs exist at a time")
                
                if max_samples:
                    logger.info(f"   Limiting to {max_samples} sample(s) for processing")
                
                start_ts = datetime.utcnow()
                stats = _steps_mod.run_download_quant_workflow(
                    metadata_path=metadata_file,
                    getfastq_params=getfastq_params,  # num_download_workers will be filtered out by the function
                    quant_params=quant_params,
                    work_dir=config.work_dir,
                    log_dir=(config.log_dir or (config.work_dir / "logs")),
                    num_workers=num_workers,  # Controls sequential (1) vs parallel (>1) mode
                    max_samples=max_samples,  # Limit number of samples to process
                )
                end_ts = datetime.utcnow()
                duration_s = max(0.0, (end_ts - start_ts).total_seconds())
                
                # Record combined getfastq+quant step
                success_rate = (stats["processed"] + stats["skipped"]) / max(1, stats["total_samples"])
                result_code = 0 if success_rate > 0.5 else 1  # Consider successful if >50% completed
                
                return_codes.append(result_code)
                manifest_records.append(
                    {
                        "step": "getfastq+quant (immediate)",
                        "return_code": result_code,
                        "stdout_bytes": 0,
                        "stderr_bytes": 0,
                        "started_utc": start_ts.isoformat() + "Z",
                        "finished_utc": end_ts.isoformat() + "Z",
                        "duration_seconds": duration_s,
                        "work_dir": str(config.work_dir),
                        "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                        "params": {
                            "getfastq": dict(getfastq_params),
                            "quant": dict(quant_params),
                            "statistics": stats,
                        },
                        "command": f"sequential_download_quant (immediate per-sample processing)",
                        "note": f"Processed {stats['processed']}, skipped {stats['skipped']}, failed {stats['failed']}/{stats['total_samples']} (immediate: download→quant→delete per sample)",
                    }
                )
        
        # Run post-steps (merge, curate, sanity, etc.)
        steps = post_steps  # Continue with post-steps in the main loop
    
    for subcommand, params in steps:
        runner = _steps_mod.STEP_RUNNERS.get(subcommand)
        if runner is None:  # pragma: no cover - defensive
            raise KeyError(f"No runner registered for step: {subcommand}")
        filtered = _sanitize_params_for_subcommand(subcommand, params)

        # Skip downstream steps when no metadata rows exist
        if not has_metadata_rows and subcommand in {
            "select",
            "getfastq",
            "integrate",
            "quant",
            "merge",
            "cstmm",
            "curate",
            "csca",
        }:
            manifest_records.append(
                {
                    "step": subcommand,
                    "return_code": 204,  # No Content
                    "stdout_bytes": 0,
                    "stderr_bytes": 0,
                    "started_utc": datetime.utcnow().isoformat() + "Z",
                    "finished_utc": datetime.utcnow().isoformat() + "Z",
                    "duration_seconds": 0.0,
                    "work_dir": str(config.work_dir),
                    "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                    "params": dict(filtered),
                    "command": " ".join(build_amalgkit_command(subcommand, filtered)),
                    "note": f"skipped: metadata has no records at {metadata_table_path}",
                }
            )
            return_codes.append(204)
            continue

        # Optional pre-filter: require tissue metadata before selection
        if subcommand == "select" and bool(config.filters.get("require_tissue", False)):
            # Infer the metadata path if not provided
            meta_path = filtered.get("metadata")
            if not meta_path or meta_path == "inferred":
                meta_path = str(config.work_dir / "metadata" / "metadata.tsv")
            meta_path_p = Path(str(meta_path))
            if meta_path_p.exists():
                out_filtered = meta_path_p.with_name("metadata.filtered.tissue.tsv")
                try:
                    rows = [
                        row
                        for row in read_delimited(meta_path_p, delimiter="\t")
                        if (row.get("tissue", "").strip() != "") and (row.get("run", "").strip() != "")
                    ]
                    write_delimited(rows, out_filtered, delimiter="\t")
                    filtered["metadata"] = str(out_filtered)
                    manifest_records.append(
                        {
                            "step": "preselect-filter",
                            "return_code": 0,
                            "stdout_bytes": 0,
                            "stderr_bytes": 0,
                            "started_utc": datetime.utcnow().isoformat() + "Z",
                            "finished_utc": datetime.utcnow().isoformat() + "Z",
                            "duration_seconds": 0.0,
                            "work_dir": str(config.work_dir),
                            "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                            "params": {"require_tissue": True, "input": str(meta_path_p), "output": str(out_filtered)},
                            "command": f"filter-metadata require_tissue=1 < {meta_path_p} > {out_filtered}",
                        }
                    )
                except Exception as exc:  # pragma: no cover - defensive
                    manifest_records.append(
                        {
                            "step": "preselect-filter",
                            "return_code": 1,
                            "stdout_bytes": 0,
                            "stderr_bytes": len(str(exc)),
                            "started_utc": datetime.utcnow().isoformat() + "Z",
                            "finished_utc": datetime.utcnow().isoformat() + "Z",
                            "duration_seconds": 0.0,
                            "work_dir": str(config.work_dir),
                            "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                            "params": {"require_tissue": True, "input": str(meta_path_p)},
                            "command": "filter-metadata require_tissue=1",
                            "note": "metadata filter failed; proceeding with original metadata",
                        }
                    )

        # Inject select.config_dir if not provided: prefer config_base produced by `amalgkit config`
        if subcommand == "select":
            cfg_dir = filtered.get("config_dir")
            if not cfg_dir or cfg_dir == "inferred":
                preferred = config.work_dir / "config_base"
                fallback = config.work_dir / "config"
                use_dir = preferred if preferred.exists() else fallback
                filtered["config_dir"] = str(use_dir)
            else:
                # Ensure config_dir is an absolute resolved path
                cfg_path = Path(str(cfg_dir)).expanduser()
                try:
                    cfg_path = cfg_path.resolve()
                except Exception:
                    # If resolution fails, try relative to work_dir
                    if not cfg_path.is_absolute():
                        cfg_path = (config.work_dir / cfg_path).resolve()
                # Verify directory exists
                if not cfg_path.exists():
                    # Try to find it relative to work_dir
                    work_relative = config.work_dir / cfg_path
                    if work_relative.exists():
                        cfg_path = work_relative.resolve()
                    else:
                        # Fallback to config_base or config
                        preferred = config.work_dir / "config_base"
                        fallback = config.work_dir / "config"
                        cfg_path = preferred if preferred.exists() else fallback
                        if not cfg_path.exists():
                            cfg_path = preferred  # Use preferred even if it doesn't exist yet
                filtered["config_dir"] = str(cfg_path)
        
        # Inject metadata path for steps that need it if not provided
        # CRITICAL: getfastq, integrate, quant, merge need row-per-sample format (NOT pivot tables)
        # Use smart fallback logic: metadata.filtered.tissue → metadata.tsv (skip pivot files)
        # select also needs metadata but can use the standard metadata.tsv format
        if subcommand == "select":
            if not filtered.get("metadata") or filtered.get("metadata") == "inferred":
                filtered["metadata"] = str(config.work_dir / "metadata" / "metadata.tsv")
        
        if subcommand in {"integrate", "getfastq", "quant", "merge"}:
            current_metadata = filtered.get("metadata")
            needs_fallback = False
            
            # Check if we need fallback: no metadata set, or existing metadata file is empty/pivot format
            if not current_metadata:
                needs_fallback = True
            else:
                # Check if current metadata file has actual data rows AND correct format (has 'run' column)
                try:
                    current_path = Path(str(current_metadata))
                    if current_path.exists():
                        rows = list(read_delimited(current_path, delimiter="\t"))
                        if len(rows) == 0:  # Only header, no data
                            needs_fallback = True
                        elif rows and 'run' not in rows[0]:  # Missing 'run' column (likely pivot table)
                            needs_fallback = True
                    else:
                        needs_fallback = True
                except Exception:
                    needs_fallback = True
            
            if needs_fallback:
                # Try to find a metadata file with actual data rows in row-per-sample format
                # SKIP pivot tables (pivot_qualified.tsv, pivot_selected.tsv) - they lack run IDs
                candidate_paths = [
                    config.work_dir / "metadata" / "metadata.filtered.tissue.tsv",
                    config.work_dir / "metadata" / "metadata.tsv",
                ]
                for candidate in candidate_paths:
                    if candidate.exists():
                        # Check if file has data rows AND 'run' column
                        try:
                            rows = list(read_delimited(candidate, delimiter="\t"))
                            if len(rows) > 0 and 'run' in rows[0]:  # Has data and correct format
                                filtered["metadata"] = str(candidate)
                                break
                        except Exception:
                            continue
        
        # Ensure fastq_dir exists before integrate step
        if subcommand == "integrate":
            fastq_dir = filtered.get("fastq_dir")
            if fastq_dir:
                # Resolve the path relative to work_dir or as absolute
                fastq_path = Path(fastq_dir)
                if not fastq_path.is_absolute():
                    # If relative, make it relative to the repository root (where we're running from)
                    fastq_path = Path.cwd() / fastq_dir
                fastq_path.mkdir(parents=True, exist_ok=True)
        
        if subcommand == "merge":
            # Ensure 'out' has a default when only out_dir provided
            if "out" not in filtered:
                default_out = config.work_dir.parent / "merged" / "merged_abundance.tsv"
                filtered["out"] = str(default_out)
        # Preflight dependency checks for steps that require extra CLIs
        dep = check_step_dependencies(subcommand)
        if not dep.ok:
            # Record skip with a conventional code (126) and continue
            manifest_records.append(
                {
                    "step": subcommand,
                    "return_code": 126,
                    "stdout_bytes": 0,
                    "stderr_bytes": 0,
                    "started_utc": datetime.utcnow().isoformat() + "Z",
                    "finished_utc": datetime.utcnow().isoformat() + "Z",
                    "duration_seconds": 0.0,
                    "work_dir": str(config.work_dir),
                    "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                    "params": dict(filtered),
                    "command": " ".join(build_amalgkit_command(subcommand, filtered)),
                    "note": f"skipped: missing dependencies: {', '.join(dep.missing)}",
                }
            )
            return_codes.append(126)
            continue

        start_ts = datetime.utcnow()
        result = runner(
            filtered,
            work_dir=config.work_dir,
            log_dir=(config.log_dir or (config.work_dir / "logs")),
            check=check,
        )
        end_ts = datetime.utcnow()
        duration_s = max(0.0, (end_ts - start_ts).total_seconds())
        return_codes.append(result.returncode)
        manifest_records.append(
            {
                "step": subcommand,
                "return_code": result.returncode,
                "stdout_bytes": len(result.stdout or ""),
                "stderr_bytes": len(result.stderr or ""),
                "started_utc": start_ts.isoformat() + "Z",
                "finished_utc": end_ts.isoformat() + "Z",
                "duration_seconds": duration_s,
                "work_dir": str(config.work_dir),
                "log_dir": str(config.log_dir or (config.work_dir / "logs")),
                "params": dict(filtered),
                "command": " ".join(build_amalgkit_command(subcommand, filtered)),
            }
        )
        if check and result.returncode != 0:
            break

    _write_manifest_records(manifest_path, manifest_records)
    _write_run_reports(config, manifest_records)
    return return_codes


def _write_run_reports(config: AmalgkitWorkflowConfig, records: list[dict[str, Any]]) -> None:
    """Emit JSON and Markdown summaries next to the manifest.

    Files:
      - amalgkit.report.json
      - amalgkit.report.md
    """
    report_json = config.work_dir / "amalgkit.report.json"
    report_md = config.work_dir / "amalgkit.report.md"
    ensure_directory(config.work_dir)

    # JSON: include config summary and records
    summary = {
        "work_dir": str(config.work_dir),
        "log_dir": str(config.log_dir or (config.work_dir / "logs")),
        "threads": config.threads,
        "species_list": list(config.species_list),
        "num_steps": len(records),
        "return_codes": [r.get("return_code", -1) for r in records],
    }
    dump_json({"summary": summary, "records": records}, report_json, indent=2)

    # Markdown: brief human-readable view
    lines: list[str] = []
    lines.append(f"# Amalgkit Run Report\n")
    lines.append(f"Work dir: `{config.work_dir}`  ")
    lines.append(f"Logs: `{config.log_dir or (config.work_dir / 'logs')}`  ")
    lines.append(f"Threads: {config.threads}  ")
    if config.species_list:
        lines.append(f"Species: {', '.join(config.species_list)}  ")
    lines.append("")
    lines.append("| Step | Code | Duration (s) |")
    lines.append("|------|------|--------------|")
    for rec in records:
        lines.append(f"| {rec['step']} | {rec['return_code']} | {rec.get('duration_seconds', 0.0):.2f} |")
    report_md.write_text("\n".join(lines), encoding="utf-8")


def load_workflow_config(config_file: str | Path) -> AmalgkitWorkflowConfig:
    """Load `AmalgkitWorkflowConfig` from a config file with env overrides.

    Paths in the config file are resolved relative to repository root.
    This allows configs to work whether the repo is on `/home/q/...` or `/media/q/ext6/...`.
    
    Expected top-level keys:
      - work_dir (str) - resolved relative to repo root
      - log_dir (str, optional) - resolved relative to repo root
      - threads (int)
      - species_list (list[str])
      - steps (mapping of step name -> params mapping)
      - auto_install_amalgkit (bool, optional)
      - genome (mapping, optional) - dest_dir resolved relative to repo root
    """
    config_path = Path(config_file).resolve()
    repo_root = _find_repo_root(config_path)
    
    raw = load_mapping_from_file(config_file)
    raw = apply_env_overrides(raw, prefix="AK")

    # Resolve work_dir relative to repo root
    work_dir_str = raw.get("work_dir", "output/amalgkit/work")
    work_dir = _resolve_path_relative_to_repo(work_dir_str, repo_root)
    
    # Resolve log_dir relative to repo root
    log_dir_val = raw.get("log_dir")
    log_dir = _resolve_path_relative_to_repo(log_dir_val, repo_root) if isinstance(log_dir_val, str) else None
    
    threads = int(raw.get("threads", 4))
    # Accept YAML list for species_list
    species_raw = raw.get("species_list", [])
    if isinstance(species_raw, list):
        species_list = [str(x) for x in species_raw]
    else:
        species_list = []
    steps_map = raw.get("steps", {}) or {}
    if not isinstance(steps_map, dict):
        steps_map = {}
    # Keep only known step names if provided
    allowed = {
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
    }
    steps_map = {k: v for k, v in steps_map.items() if str(k) in allowed and isinstance(v, dict)}
    
    # Resolve paths in step configs relative to repo root
    for step_name, step_params in steps_map.items():
        if isinstance(step_params, dict):
            # Resolve out_dir if present
            if "out_dir" in step_params:
                step_params["out_dir"] = str(_resolve_path_relative_to_repo(step_params["out_dir"], repo_root))
            # Resolve out if present (for merge step)
            if "out" in step_params:
                step_params["out"] = str(_resolve_path_relative_to_repo(step_params["out"], repo_root))

    auto_install_amalgkit = bool(raw.get("auto_install_amalgkit", False))
    genome_cfg = raw.get("genome") if isinstance(raw.get("genome"), dict) else None
    
    # Resolve genome dest_dir relative to repo root
    if genome_cfg and "dest_dir" in genome_cfg:
        genome_cfg["dest_dir"] = str(_resolve_path_relative_to_repo(genome_cfg["dest_dir"], repo_root))
    
    filters = raw.get("filters") if isinstance(raw.get("filters"), dict) else {}

    return AmalgkitWorkflowConfig(
        work_dir=work_dir,
        log_dir=log_dir,
        threads=threads,
        species_list=species_list,
        per_step={str(k): dict(v) for k, v in steps_map.items()},
        auto_install_amalgkit=auto_install_amalgkit,
        genome=genome_cfg,
        filters=filters,
    )


__all__ = [
    "AmalgkitWorkflowConfig",
    "plan_workflow",
    "execute_workflow",
    "load_workflow_config",
]
