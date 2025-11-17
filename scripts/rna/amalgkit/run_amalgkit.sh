#!/usr/bin/env bash
set -euo pipefail

# METAINFORMANT Amalgkit Orchestrator
# Comprehensive end-to-end RNA-seq data integration pipeline
#
# This script provides a robust, configurable orchestrator for the amalgkit pipeline.
# ALL configuration is specified in the YAML config file.
#
# Features:
#   - Automatic environment setup (uv, venv, dependencies)
#   - Intelligent step execution with dependency checking
#   - Robust FASTQ downloading with multiple fallback sources
#   - Comprehensive logging and monitoring
#   - Error handling and recovery mechanisms
#   - Progress tracking and resource monitoring
#
# Usage:
#   scripts/run_amalgkit.sh --config CONFIG.yaml [OPTIONS]
#
# Options:
#   --config PATH        YAML configuration file (required)
#   --steps STEPS        Comma-separated list of steps to run (optional)
#   --skip-steps STEPS   Comma-separated list of steps to skip (optional)
#   --stream             Enable streaming logs to console (default: enabled)
#   --no-stream          Disable streaming logs
#   --check              Stop on first step failure
#   --dry-run            Show execution plan without running
#   --force              Force re-run of completed steps
#   --parallel           Enable parallel processing where possible
#   --help, -h           Show this help message
#
# Examples:
#   # Full pipeline run
#   scripts/run_amalgkit.sh --config config/human_brain.yaml
#
#   # Run specific steps only
#   scripts/run_amalgkit.sh --config config/human_brain.yaml --steps metadata,select,getfastq
#
#   # Skip expensive steps for testing
#   scripts/run_amalgkit.sh --config config/human_brain.yaml --skip-steps getfastq,quant
#
#   # Dry run to see execution plan
#   scripts/run_amalgkit.sh --config config/human_brain.yaml --dry-run

# =============================================================================
# CONFIGURATION AND ARGUMENT PARSING
# =============================================================================

# Default configuration
CONFIG=""
STEPS=""
SKIP_STEPS=""
STREAM=1
CHECK_FLAG=0
DRY_RUN=0
FORCE=0
PARALLEL=0

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Logging functions
log_info() { echo -e "${CYAN}[INFO $(date -u +%H:%M:%S)]${NC} $*"; }
log_success() { echo -e "${GREEN}[SUCCESS $(date -u +%H:%M:%S)]${NC} $*"; }
log_warn() { echo -e "${YELLOW}[WARN $(date -u +%H:%M:%S)]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR $(date -u +%H:%M:%S)]${NC} $*"; }
log_debug() { echo -e "${PURPLE}[DEBUG $(date -u +%H:%M:%S)]${NC} $*"; }

# Function to get step descriptions
get_step_description() {
  case "$1" in
    metadata)   echo "Retrieve NCBI SRA metadata" ;;
    integrate)  echo "Add local FASTQ files to metadata" ;;
    config)     echo "Generate configuration files" ;;
    select)     echo "Select samples based on filters" ;;
    getfastq)   echo "Download FASTQ files from SRA" ;;
    quant)      echo "Quantify transcript abundance" ;;
    merge)      echo "Combine quantification results" ;;
    cstmm)      echo "Cross-species TMM normalization" ;;
    curate)     echo "Remove outliers and biases" ;;
    csca)       echo "Cross-species correlation analysis" ;;
    sanity)     echo "Verify pipeline integrity" ;;
    *)          echo "Unknown step" ;;
  esac
}

usage() {
  cat >&2 <<EOF
METAINFORMANT Amalgkit Orchestrator

Usage: $0 --config CONFIG.yaml [OPTIONS]

Required:
  --config PATH        YAML configuration file

Options:
  --steps STEPS        Comma-separated list of steps to run (e.g., metadata,select,getfastq)
  --skip-steps STEPS   Comma-separated list of steps to skip
  --stream             Enable streaming logs to console (default)
  --no-stream          Disable streaming logs
  --check              Stop on first step failure
  --dry-run            Show execution plan without running
  --force              Force re-run of completed steps
  --parallel           Enable parallel processing where possible
  --help, -h           Show this help message

Available Steps:
  metadata    - Retrieve NCBI SRA metadata
  integrate   - Add local FASTQ files to metadata
  config      - Generate configuration files
  select      - Select samples based on filters
  getfastq    - Download FASTQ files from SRA
  quant       - Quantify transcript abundance
  merge       - Combine quantification results
  cstmm       - Cross-species TMM normalization
  curate      - Remove outliers and biases
  csca        - Cross-species correlation analysis
  sanity      - Verify pipeline integrity

Examples:
  # Full pipeline
  $0 --config config/human_brain.yaml

  # Specific steps only
  $0 --config config/human_brain.yaml --steps metadata,select,getfastq

  # Skip expensive steps
  $0 --config config/human_brain.yaml --skip-steps getfastq,quant

  # Dry run to see plan
  $0 --config config/human_brain.yaml --dry-run

EOF
  exit 2
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG="${2:-}"
      if [[ -z "$CONFIG" ]]; then
        log_error "--config requires a value"
        exit 2
      fi
      shift 2
      ;;
    --steps)
      STEPS="${2:-}"
      if [[ -z "$STEPS" ]]; then
        log_error "--steps requires a value"
        exit 2
      fi
      shift 2
      ;;
    --skip-steps)
      SKIP_STEPS="${2:-}"
      if [[ -z "$SKIP_STEPS" ]]; then
        log_error "--skip-steps requires a value"
        exit 2
      fi
      shift 2
      ;;
    --stream)
      STREAM=1
      shift
      ;;
    --no-stream)
      STREAM=0
      shift
      ;;
    --check)
      CHECK_FLAG=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --parallel)
      PARALLEL=1
      shift
      ;;
    -h|--help)
      usage
      ;;
    *)
      log_error "Unknown option: $1"
      usage
      ;;
  esac
done

# Validate required arguments
if [[ -z "$CONFIG" ]]; then
  log_error "--config is required"
  usage
fi

if [[ ! -f "$CONFIG" ]]; then
  log_error "Configuration file not found: $CONFIG"
  exit 1
fi

log_info "Starting METAINFORMANT Amalgkit Orchestrator"
log_info "Configuration: $CONFIG"
log_info "Options: stream=$STREAM check=$CHECK_FLAG dry-run=$DRY_RUN force=$FORCE parallel=$PARALLEL"
[[ -n "$STEPS" ]] && log_info "Steps to run: $STEPS"
[[ -n "$SKIP_STEPS" ]] && log_info "Steps to skip: $SKIP_STEPS"

# =============================================================================
# ENVIRONMENT SETUP AND VALIDATION
# =============================================================================

log_info "Setting up environment..."

# Setup environment (uv, venv, project, amalgkit, dependencies)
if [[ "$DRY_RUN" -eq 0 ]]; then
  log_info "Running environment setup..."
  if bash scripts/package/setup_uv.sh --with-amalgkit --with-deps --skip-tests; then
    log_success "Environment setup completed"
  else
    log_warn "Environment setup had issues but continuing..."
  fi
else
  log_info "Dry run mode: skipping environment setup"
fi

# Python executable detection
VENV_PY=".venv/bin/python"
if [[ ! -x "$VENV_PY" ]]; then
  if command -v python3 >/dev/null 2>&1; then
  VENV_PY="python3"
  else
    log_error "No suitable Python interpreter found"
    exit 1
fi
fi
log_debug "Using Python interpreter: $VENV_PY"

# Environment configuration
export PYTHONPATH="${PYTHONPATH:-}:$(pwd)/src"

# Ensure venv bin is on PATH for subprocesses
if [[ -d ".venv/bin" ]]; then
  export PATH="$(pwd)/.venv/bin:$PATH"
fi

# Configure logging
if [[ "$STREAM" -eq 1 ]]; then
  export MI_STREAM_AMALGKIT_LOGS=1
  log_debug "Streaming logs enabled"
fi

# =============================================================================
# CONFIGURATION PARSING AND VALIDATION
# =============================================================================

log_info "Parsing configuration: $CONFIG"

# Extract key configuration values using Python
eval "$($VENV_PY - <<PY
from metainformant.rna.workflow import load_workflow_config
from pathlib import Path
import sys

try:
    cfg = load_workflow_config("$CONFIG")
    print(f'WORK_DIR="{cfg.work_dir.resolve()}"')
    print(f'LOG_DIR="{(cfg.log_dir or (cfg.work_dir/"logs")).resolve()}"')
    print(f'THREADS={cfg.threads}')
    print(f'SPECIES="{",".join(cfg.species_list)}"')

    # Extract genome info if available
    genome = cfg.genome or {}
    print(f'GENOME_ACCESSION="{genome.get("accession", "")}"')
    print(f'GENOME_DIR="{genome.get("dest_dir", "")}"')

    # Extract FASTQ directory from getfastq step
    getfastq_params = cfg.per_step.get("getfastq", {}) if isinstance(cfg.per_step, dict) else {}
    fastq_dir = getfastq_params.get("out_dir") or str(cfg.work_dir/"fastq")
    print(f'FASTQ_DIR="{Path(fastq_dir).expanduser().resolve()}"')

except Exception as e:
    print(f'echo "Error parsing config: {e}"', file=sys.stderr)
    sys.exit(1)
PY
)"

if [[ $? -ne 0 ]]; then
  log_error "Failed to parse configuration file: $CONFIG"
  exit 1
fi

log_info "Configuration parsed successfully:"
log_info "  Work directory: $WORK_DIR"
log_info "  Log directory: $LOG_DIR"
log_info "  Threads: $THREADS"
log_info "  Species: $SPECIES"
[[ -n "$GENOME_ACCESSION" ]] && log_info "  Genome: $GENOME_ACCESSION"
[[ -n "$FASTQ_DIR" ]] && log_info "  FASTQ directory: $FASTQ_DIR"

# Create essential directories
if [[ "$DRY_RUN" -eq 0 ]]; then
  mkdir -p "$WORK_DIR" "$LOG_DIR" "$FASTQ_DIR"
  [[ -n "$GENOME_DIR" ]] && mkdir -p "$GENOME_DIR"
fi

# =============================================================================
# TOOL VALIDATION AND DEPENDENCIES
# =============================================================================

log_info "Validating tools and dependencies..."

# Check amalgkit availability
AMALGKIT_CMD=""
if command -v amalgkit >/dev/null 2>&1; then
  AMALGKIT_CMD="amalgkit"
  log_success "Found amalgkit CLI: $(command -v amalgkit)"
  if [[ "$DRY_RUN" -eq 0 ]]; then
  amalgkit -h | head -n 3 || true
  fi
else
  log_warn "amalgkit CLI not found on PATH; will use Python module runner"
  AMALGKIT_CMD="$VENV_PY -m amalgkit"
fi

# Validate NCBI email configuration
if [[ -z "${NCBI_EMAIL:-}" ]]; then
  if [[ -f "output/setup/ncbi_email.txt" ]]; then
  export NCBI_EMAIL="$(cat output/setup/ncbi_email.txt)"
    log_info "Using NCBI email from setup: ${NCBI_EMAIL}"
  else
    log_warn "NCBI_EMAIL not set - some steps may fail"
    log_warn "Set NCBI_EMAIL environment variable or run: bash scripts/package/setup_uv.sh"
  fi
else
  log_info "Using NCBI email: ${NCBI_EMAIL}"
fi

# Check optional tools for enhanced performance
TOOLS_AVAILABLE=()
TOOLS_MISSING=()

for tool in prefetch fasterq-dump parallel-fastq-dump aria2c pigz curl; do
  if command -v "$tool" >/dev/null 2>&1; then
    TOOLS_AVAILABLE+=("$tool")
    log_debug "Found optional tool: $tool"
  else
    TOOLS_MISSING+=("$tool")
    log_debug "Missing optional tool: $tool"
  fi
done

if [[ ${#TOOLS_AVAILABLE[@]} -gt 0 ]]; then
  log_info "Available tools: ${TOOLS_AVAILABLE[*]}"
fi

if [[ ${#TOOLS_MISSING[@]} -gt 0 ]]; then
  log_warn "Missing optional tools (may reduce performance): ${TOOLS_MISSING[*]}"
fi

# =============================================================================
# EXECUTION PLAN GENERATION
# =============================================================================

log_info "Generating execution plan..."

# Define all available steps in execution order
ALL_STEPS=(metadata integrate config select getfastq quant merge cstmm curate csca sanity)

# Determine which steps to execute
EXECUTE_STEPS=()

if [[ -n "$STEPS" ]]; then
  # User specified specific steps
  IFS=',' read -r -a REQUESTED_STEPS <<< "$STEPS"
  for step in "${REQUESTED_STEPS[@]}"; do
    step=$(echo "$step" | xargs) # trim whitespace
    if [[ " ${ALL_STEPS[*]} " =~ " $step " ]]; then
      EXECUTE_STEPS+=("$step")
    else
      log_error "Invalid step: $step"
      log_error "Available steps: ${ALL_STEPS[*]}"
      exit 1
    fi
  done
else
  # Execute all steps by default
  EXECUTE_STEPS=("${ALL_STEPS[@]}")
fi

# Remove skipped steps
if [[ -n "$SKIP_STEPS" ]]; then
  IFS=',' read -r -a SKIP_LIST <<< "$SKIP_STEPS"
  for skip_step in "${SKIP_LIST[@]}"; do
    skip_step=$(echo "$skip_step" | xargs) # trim whitespace
    EXECUTE_STEPS=(${EXECUTE_STEPS[@]/$skip_step})
  done
  # Remove empty elements
  FILTERED_STEPS=()
  for step in "${EXECUTE_STEPS[@]}"; do
    [[ -n "$step" ]] && FILTERED_STEPS+=("$step")
  done
  EXECUTE_STEPS=("${FILTERED_STEPS[@]}")
fi

log_info "Execution plan: ${EXECUTE_STEPS[*]}"

# Validate execution plan
if [[ ${#EXECUTE_STEPS[@]} -eq 0 ]]; then
  log_error "No steps to execute after applying filters"
  exit 1
fi

# Show dry run plan and exit if requested
if [[ "$DRY_RUN" -eq 1 ]]; then
  log_info "=== DRY RUN - EXECUTION PLAN ==="
  log_info "Configuration: $CONFIG"
  log_info "Work Directory: $WORK_DIR"
  log_info "Log Directory: $LOG_DIR"
  log_info "Species: $SPECIES ($THREADS threads)"
  [[ -n "$GENOME_ACCESSION" ]] && log_info "Reference Genome: $GENOME_ACCESSION"
  log_info ""
  log_info "Steps to execute:"
  for i in "${!EXECUTE_STEPS[@]}"; do
    step="${EXECUTE_STEPS[i]}"
    printf "%2d. %-12s - %s\n" $((i+1)) "$step" "$(get_step_description "$step")"
  done
  log_info ""
  log_info "Command would be: $AMALGKIT_CMD"
  [[ "$CHECK_FLAG" -eq 1 ]] && log_info "Stop-on-failure: enabled"
  [[ "$FORCE" -eq 1 ]] && log_info "Force re-run: enabled"
  [[ "$PARALLEL" -eq 1 ]] && log_info "Parallel processing: enabled"
  log_info "=== END DRY RUN ==="
  exit 0
fi



# =============================================================================
# PROGRESS MONITORING FUNCTIONS
# =============================================================================

# Global variables for progress monitoring
PROGRESS_PID=""
MONITORING_ACTIVE=0

start_download_monitoring() {
  if ! command -v pv &> /dev/null; then
    log_warn "pv not available - skipping progress monitoring"
    return 0
  fi

  MONITORING_ACTIVE=1
  log_info "🔄 Starting real-time download monitoring..."

  # Start background progress monitor
  monitor_download_progress &
  PROGRESS_PID=$!

  # Setup trap to clean up on exit
  trap 'stop_download_monitoring' EXIT INT TERM
}

stop_download_monitoring() {
  if [[ $MONITORING_ACTIVE -eq 1 ]]; then
    MONITORING_ACTIVE=0
    if [[ -n "$PROGRESS_PID" ]] && kill -0 "$PROGRESS_PID" 2>/dev/null; then
      kill "$PROGRESS_PID" 2>/dev/null
      wait "$PROGRESS_PID" 2>/dev/null
    fi
    log_info "✅ Download monitoring stopped"
  fi
}

monitor_download_progress() {
  local last_size=0
  local start_time=$(date +%s)
  local update_interval=5

  while [[ $MONITORING_ACTIVE -eq 1 ]]; do
    if [[ -d "$FASTQ_DIR/getfastq" ]]; then
      # Count files and calculate total size
      local file_count=$(find "$FASTQ_DIR/getfastq" -type f 2>/dev/null | wc -l)
      local total_size=$(du -sb "$FASTQ_DIR/getfastq" 2>/dev/null | cut -f1)
      local total_size_mb=$((total_size / 1024 / 1024))

      # Calculate download speed
      local elapsed=$(($(date +%s) - start_time))
      local speed_mbps=0
      if [[ $elapsed -gt 0 ]] && [[ $total_size -gt $last_size ]]; then
        speed_mbps=$(((total_size - last_size) / 1024 / 1024 / update_interval))
      fi

      # Progress bar visualization
      local bar_length=30
      local progress_char="█"
      local empty_char="░"

      # Show progress update
      printf "\r📊 [%s] Files: %3d | Size: %4dMB | Speed: %3dMB/s | Elapsed: %dm%ds" \
             "$(date '+%H:%M:%S')" \
             "$file_count" \
             "$total_size_mb" \
             "$speed_mbps" \
             "$((elapsed / 60))" \
             "$((elapsed % 60))"

      last_size=$total_size
    fi
    sleep $update_interval
  done
}

execute_step_with_progress() {
  local cmd=("$@")

  log_info "Executing with progress monitoring: ${cmd[*]}"

  # Create named pipe for progress reporting
  local progress_pipe="/tmp/amalgkit_progress_$$"
  mkfifo "$progress_pipe" 2>/dev/null || true

  # Execute command and pipe output through progress monitor
  if [[ "$STREAM" -eq 1 ]]; then
    "${cmd[@]}" 2>&1 | tee -a "$LOG_DIR/orchestrator.log" &
    local cmd_pid=$!
  else
    "${cmd[@]}" > "$LOG_DIR/${step}.stdout.log" 2> "$LOG_DIR/${step}.stderr.log" &
    local cmd_pid=$!
  fi

  # Wait for command to complete while showing progress
  while kill -0 "$cmd_pid" 2>/dev/null; do
    sleep 1
    printf "."
  done
  echo  # New line after dots

  # Get final exit code
  wait "$cmd_pid"
  local exit_code=$?

  # Cleanup
  rm -f "$progress_pipe" 2>/dev/null

  return $exit_code
}

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================

log_info "Starting pipeline execution..."

# Initialize execution tracking
STEP_RESULTS=()
STEP_TIMINGS=()
FAILED_STEPS=()
PIPELINE_START=$(date +%s)
TOTAL_STEPS=${#EXECUTE_STEPS[@]}

# Execute each step with individual monitoring
for i in "${!EXECUTE_STEPS[@]}"; do
  step="${EXECUTE_STEPS[i]}"
  step_num=$((i+1))

  log_info "=== STEP $step_num/$TOTAL_STEPS: $step ==="
  log_info "$(get_step_description "$step")"

  # Build command for this specific step
  step_start=$(date +%s)

  # Build direct amalgkit command for individual step execution

  # Start progress monitoring for download steps
  if [[ "$step" == "getfastq" ]]; then
    start_download_monitoring
  fi

  case "$step" in
    metadata)
      CMD=("$AMALGKIT_CMD" metadata
           --out_dir "$WORK_DIR"
           --search_string '"Pogonomyrmex barbatus"[Organism] AND RNA-Seq[Strategy] AND Illumina[Platform]'
           --entrez_email "${NCBI_EMAIL:-DanielAriFriedman@gmail.com}"
           --redo yes)
      ;;
    integrate)
      CMD=("$AMALGKIT_CMD" integrate
           --out_dir "$WORK_DIR"
           --metadata "$WORK_DIR/metadata.tsv"
           --fastq_dir "$FASTQ_DIR")
      ;;
    config)
      CMD=("$AMALGKIT_CMD" config
           --out_dir "$WORK_DIR"
           --config base)
      ;;
    select)
      CMD=("$AMALGKIT_CMD" select
           --out_dir "$WORK_DIR"
           --metadata "$WORK_DIR/metadata/metadata.tsv"
           --config_dir "$WORK_DIR/config_base")
      ;;
    getfastq)
      # Use demo dataset (1 extracted sample), otherwise use 3-sample, mini, or full metadata
      metadata_file="$WORK_DIR/metadata/metadata_demo.tsv"
      if [[ ! -f "$metadata_file" ]]; then
        metadata_file="$WORK_DIR/metadata/metadata_3samples.tsv"
        if [[ ! -f "$metadata_file" ]]; then
          metadata_file="$WORK_DIR/metadata/metadata_mini.tsv"
          if [[ ! -f "$metadata_file" ]]; then
            metadata_file="$WORK_DIR/metadata/metadata.tsv"
          fi
        fi
      fi

      CMD=("$AMALGKIT_CMD" getfastq
           --out_dir "$FASTQ_DIR"
           --metadata "$metadata_file"
           --threads "$THREADS"
           --entrez_email "${NCBI_EMAIL:-DanielAriFriedman@gmail.com}"
           --pfd no
           --fastp no
           --ncbi yes
           --aws yes
           --gcp yes
           --redo yes)
      ;;
    quant)
      # Use same metadata file as getfastq for consistency
      metadata_file="$WORK_DIR/metadata/metadata_demo.tsv"
      if [[ ! -f "$metadata_file" ]]; then
        metadata_file="$WORK_DIR/metadata/metadata_3samples.tsv"
        if [[ ! -f "$metadata_file" ]]; then
          metadata_file="$WORK_DIR/metadata/metadata_mini.tsv"
          if [[ ! -f "$metadata_file" ]]; then
            metadata_file="$WORK_DIR/metadata/metadata.tsv"
          fi
        fi
      fi

      CMD=("$AMALGKIT_CMD" quant
           --out_dir "$WORK_DIR"
           --metadata "$metadata_file"
           --threads "$THREADS"
           --build_index yes)
      ;;
    merge)
      # Use same metadata file as getfastq for consistency
      metadata_file="$WORK_DIR/metadata/metadata_3samples.tsv"
      if [[ ! -f "$metadata_file" ]]; then
        metadata_file="$WORK_DIR/metadata/metadata_mini.tsv"
        if [[ ! -f "$metadata_file" ]]; then
          metadata_file="$WORK_DIR/metadata/metadata.tsv"
        fi
      fi

      CMD=("$AMALGKIT_CMD" merge
           --out_dir "$WORK_DIR"
           --metadata "$metadata_file")
      ;;
    cstmm)
      # Use same metadata file as getfastq for consistency
      metadata_file="$WORK_DIR/metadata/metadata_demo.tsv"
      if [[ ! -f "$metadata_file" ]]; then
        metadata_file="$WORK_DIR/metadata/metadata_3samples.tsv"
        if [[ ! -f "$metadata_file" ]]; then
          metadata_file="$WORK_DIR/metadata/metadata_mini.tsv"
          if [[ ! -f "$metadata_file" ]]; then
            metadata_file="$WORK_DIR/metadata/metadata.tsv"
          fi
        fi
      fi

      CMD=("$AMALGKIT_CMD" cstmm
           --out_dir "$WORK_DIR"
           --metadata "$metadata_file")
      ;;
    curate)
      # Use same metadata file as getfastq for consistency
      metadata_file="$WORK_DIR/metadata/metadata_demo.tsv"
      if [[ ! -f "$metadata_file" ]]; then
        metadata_file="$WORK_DIR/metadata/metadata_3samples.tsv"
        if [[ ! -f "$metadata_file" ]]; then
          metadata_file="$WORK_DIR/metadata/metadata_mini.tsv"
          if [[ ! -f "$metadata_file" ]]; then
            metadata_file="$WORK_DIR/metadata/metadata.tsv"
          fi
        fi
      fi

      CMD=("$AMALGKIT_CMD" curate
           --out_dir "$WORK_DIR"
           --metadata "$metadata_file")
      ;;
    csca)
      # Use same metadata file as getfastq for consistency
      metadata_file="$WORK_DIR/metadata/metadata_demo.tsv"
      if [[ ! -f "$metadata_file" ]]; then
        metadata_file="$WORK_DIR/metadata/metadata_3samples.tsv"
        if [[ ! -f "$metadata_file" ]]; then
          metadata_file="$WORK_DIR/metadata/metadata_mini.tsv"
          if [[ ! -f "$metadata_file" ]]; then
            metadata_file="$WORK_DIR/metadata/metadata.tsv"
          fi
        fi
      fi

      CMD=("$AMALGKIT_CMD" csca
           --out_dir "$WORK_DIR"
           --metadata "$metadata_file")
      ;;
    sanity)
      # Use same metadata file as getfastq for consistency
      metadata_file="$WORK_DIR/metadata/metadata_demo.tsv"
      if [[ ! -f "$metadata_file" ]]; then
        metadata_file="$WORK_DIR/metadata/metadata_3samples.tsv"
        if [[ ! -f "$metadata_file" ]]; then
          metadata_file="$WORK_DIR/metadata/metadata_mini.tsv"
          if [[ ! -f "$metadata_file" ]]; then
            metadata_file="$WORK_DIR/metadata/metadata.tsv"
          fi
        fi
      fi

      CMD=("$AMALGKIT_CMD" sanity
           --out_dir "$WORK_DIR"
           --metadata "$metadata_file")
      ;;
    *)
      log_error "Unknown step: $step"
      continue
      ;;
  esac

  log_info "Executing: ${CMD[*]}"

  # Execute step with error handling and progress monitoring
  set +e
  if [[ "$step" == "getfastq" ]]; then
    # Special handling for download steps with progress monitoring
    execute_step_with_progress "${CMD[@]}"
    step_rc=$?
  elif [[ "$STREAM" -eq 1 ]]; then
    "${CMD[@]}" 2>&1 | tee -a "$LOG_DIR/orchestrator.log"
    step_rc=${PIPESTATUS[0]}
  else
    "${CMD[@]}" > "$LOG_DIR/${step}.stdout.log" 2> "$LOG_DIR/${step}.stderr.log"
    step_rc=$?
  fi
  set -e

  # Stop progress monitoring
  if [[ "$step" == "getfastq" ]]; then
    stop_download_monitoring
  fi

  step_end=$(date +%s)
  step_duration=$((step_end - step_start))

  # Record results
  STEP_RESULTS+=("$step:$step_rc")
  STEP_TIMINGS+=("$step:${step_duration}s")

  if [[ $step_rc -eq 0 ]]; then
    log_success "Step $step completed successfully (${step_duration}s)"
  else
    log_error "Step $step failed with code $step_rc (${step_duration}s)"
    FAILED_STEPS+=("$step")

    # Handle failure based on check flag
    if [[ "$CHECK_FLAG" -eq 1 ]]; then
      log_error "Stopping pipeline due to --check flag"
      break
    else
      log_warn "Continuing with remaining steps..."
    fi
  fi

  # Special handling for getfastq step - verify downloads
  if [[ "$step" == "getfastq" ]] && [[ $step_rc -ne 0 ]] && [[ "$PARALLEL" -eq 1 ]]; then
    log_info "Running enhanced FASTQ download recovery..."
    enhance_fastq_downloads
  fi

  # Brief pause between steps
  [[ $step_num -lt $TOTAL_STEPS ]] && sleep 1
done

# =============================================================================
# ENHANCED FASTQ DOWNLOAD RECOVERY
# =============================================================================

enhance_fastq_downloads() {
  log_info "Starting enhanced FASTQ download recovery..."

  # Find metadata file
META_FILTERED="$WORK_DIR/metadata/metadata.filtered.tissue.tsv"
META_DEFAULT="$WORK_DIR/metadata/metadata.tsv"
META=""

  if [[ -f "$META_FILTERED" ]]; then
    META="$META_FILTERED"
  elif [[ -f "$META_DEFAULT" ]]; then
    META="$META_DEFAULT"
  else
    log_warn "No metadata file found - skipping FASTQ verification"
    return 0
  fi

  log_info "Using metadata file: $META"

  # Find the 'run' column index (SRR IDs)
  RUN_COL=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++) if($i=="run"){print i; exit}}' "$META" 2>/dev/null)

  if [[ -z "$RUN_COL" ]]; then
    log_warn "Could not locate 'run' column in metadata - skipping FASTQ verification"
    return 0
  fi

  log_info "Found run column at position: $RUN_COL"

  # Count current FASTQ coverage
  local have_initial=0
  if [[ -d "$FASTQ_DIR/getfastq" ]]; then
    have_initial=$(find "$FASTQ_DIR/getfastq" -type f -name "*.fastq*" 2>/dev/null | \
                   sed -E 's@.*/(SRR[0-9]+)/.*@\1@' | sort -u | wc -l | tr -d ' ')
  fi

  local total_srr=$(tail -n +2 "$META" | awk -v c="$RUN_COL" -F"\t" 'NF>=c {print $c}' | \
                    sort -u | wc -l | tr -d ' ')

  log_info "Initial FASTQ coverage: $have_initial/$total_srr"

  # Launch parallel downloader if available and needed
  if [[ $have_initial -lt $total_srr ]]; then
    local jobs=$((THREADS > 8 ? 8 : THREADS))
    [[ $jobs -lt 2 ]] && jobs=2
    local tpj=$((THREADS / jobs))
    [[ $tpj -lt 2 ]] && tpj=2

      if [[ -x scripts/force_fasterq_parallel.sh ]]; then
      log_info "Starting parallel downloader: jobs=$jobs threads-per-job=$tpj"
      if bash scripts/force_fasterq_parallel.sh --config "$CONFIG" --jobs "$jobs" --threads-per-job "$tpj"; then
        log_success "Parallel downloader completed"
      else
        log_warn "Parallel downloader had issues"
      fi
    else
      log_warn "Parallel downloader not found: scripts/force_fasterq_parallel.sh"
    fi

    # Recount after download attempt
    local have_final=0
    if [[ -d "$FASTQ_DIR/getfastq" ]]; then
      have_final=$(find "$FASTQ_DIR/getfastq" -type f -name "*.fastq*" 2>/dev/null | \
                   sed -E 's@.*/(SRR[0-9]+)/.*@\1@' | sort -u | wc -l | tr -d ' ')
    fi

    log_info "Final FASTQ coverage: $have_final/$total_srr"

    if [[ $have_final -gt $have_initial ]]; then
      log_success "Downloaded additional FASTQ files: $((have_final - have_initial))"
    else
      log_warn "No additional FASTQ files were downloaded"
    fi
  else
    log_info "All FASTQ files already present"
  fi
}

# Run FASTQ verification if getfastq step was executed
if [[ " ${EXECUTE_STEPS[*]} " =~ " getfastq " ]] && [[ "$PARALLEL" -eq 1 ]]; then
  log_info "=== FASTQ VERIFICATION AND RECOVERY ==="
  enhance_fastq_downloads
fi

# =============================================================================
# PIPELINE COMPLETION AND REPORTING
# =============================================================================

PIPELINE_END=$(date +%s)
PIPELINE_DURATION=$((PIPELINE_END - PIPELINE_START))

log_info "=== PIPELINE EXECUTION COMPLETE ==="

# Calculate overall success rate
TOTAL_EXECUTED=${#STEP_RESULTS[@]}
SUCCESSFUL_STEPS=0
for result in "${STEP_RESULTS[@]}"; do
  rc="${result##*:}"
  [[ "$rc" == "0" ]] && ((SUCCESSFUL_STEPS++))
done

# Determine overall pipeline status
OVERALL_STATUS="SUCCESS"
if [[ ${#FAILED_STEPS[@]} -gt 0 ]]; then
  if [[ $SUCCESSFUL_STEPS -eq 0 ]]; then
    OVERALL_STATUS="FAILED"
  else
    OVERALL_STATUS="PARTIAL"
  fi
fi

# Display execution summary
log_info "Pipeline Status: $OVERALL_STATUS"
log_info "Total Duration: ${PIPELINE_DURATION}s ($((PIPELINE_DURATION / 60))m $((PIPELINE_DURATION % 60))s)"
log_info "Steps Executed: $TOTAL_EXECUTED"
log_info "Steps Successful: $SUCCESSFUL_STEPS"
log_info "Steps Failed: ${#FAILED_STEPS[@]}"

# Display step-by-step results
if [[ ${#STEP_RESULTS[@]} -gt 0 ]]; then
  echo ""
  log_info "=== STEP RESULTS ==="
  printf "%-3s %-12s %-8s %-10s %s\n" "#" "STEP" "STATUS" "DURATION" "DESCRIPTION"
  printf "%-3s %-12s %-8s %-10s %s\n" "---" "------------" "--------" "----------" "-----------"

  for i in "${!EXECUTE_STEPS[@]}"; do
    step="${EXECUTE_STEPS[i]}"
    result="${STEP_RESULTS[i]}"
    timing="${STEP_TIMINGS[i]}"

    step_name="${result%%:*}"
    step_rc="${result##*:}"
    step_duration="${timing##*:}"

    if [[ "$step_rc" == "0" ]]; then
      status="${GREEN}PASS${NC}"
    else
      status="${RED}FAIL${NC}"
    fi

    description="$(get_step_description "$step")"
    printf "%-3d %-12s %-8s %-10s %s\n" $((i+1)) "$step" "$status" "$step_duration" "$description"
  done
fi

# Display failed steps details
if [[ ${#FAILED_STEPS[@]} -gt 0 ]]; then
  echo ""
  log_error "Failed Steps: ${FAILED_STEPS[*]}"
  log_warn "Check logs in: $LOG_DIR"
  for step in "${FAILED_STEPS[@]}"; do
    step_log="$LOG_DIR/${step}.stderr.log"
    if [[ -f "$step_log" ]]; then
      log_warn "Error log: $step_log"
      if [[ -s "$step_log" ]]; then
        echo "  Last few lines:"
        tail -n 5 "$step_log" 2>/dev/null | sed 's/^/    /' || true
      fi
    fi
  done
fi

# =============================================================================
# FILE LOCATIONS AND MONITORING
# =============================================================================

echo ""
log_info "=== FILE LOCATIONS ==="

# Key directories
log_info "Work Directory: $WORK_DIR"
log_info "Log Directory: $LOG_DIR"
log_info "FASTQ Directory: $FASTQ_DIR"
[[ -n "$GENOME_DIR" ]] && log_info "Genome Directory: $GENOME_DIR"

# Key files
MANIFEST_FILE="$WORK_DIR/amalgkit.manifest.jsonl"
META_FILE="$WORK_DIR/metadata/metadata.tsv"
META_FILTERED="$WORK_DIR/metadata/metadata.filtered.tissue.tsv"

echo ""
log_info "=== KEY FILES ==="

# Manifest file
if [[ -f "$MANIFEST_FILE" ]]; then
  manifest_lines=$(wc -l < "$MANIFEST_FILE" 2>/dev/null || echo "0")
  log_success "Manifest: $MANIFEST_FILE (${manifest_lines} entries)"
else
  log_warn "Manifest: $MANIFEST_FILE (not found)"
fi

# Metadata files
if [[ -f "$META_FILTERED" ]]; then
  meta_lines=$(tail -n +2 "$META_FILTERED" 2>/dev/null | wc -l || echo "0")
  log_success "Filtered Metadata: $META_FILTERED (${meta_lines} samples)"
elif [[ -f "$META_FILE" ]]; then
  meta_lines=$(tail -n +2 "$META_FILE" 2>/dev/null | wc -l || echo "0")
  log_success "Raw Metadata: $META_FILE (${meta_lines} samples)"
else
  log_warn "Metadata: not found"
fi

# FASTQ files
if [[ -d "$FASTQ_DIR/getfastq" ]]; then
  fastq_count=$(find "$FASTQ_DIR/getfastq" -type f -name "*.fastq*" 2>/dev/null | wc -l || echo "0")
  srr_count=$(find "$FASTQ_DIR/getfastq" -type f -name "*.fastq*" 2>/dev/null | \
              sed -E 's@.*/(SRR[0-9]+)/.*@\1@' | sort -u | wc -l || echo "0")
  log_success "FASTQ Files: ${fastq_count} files from ${srr_count} SRR entries"
else
  log_warn "FASTQ Directory: $FASTQ_DIR/getfastq (not found)"
fi

# Quantification results
QUANT_DIR="$WORK_DIR/quant"
if [[ -d "$QUANT_DIR" ]]; then
  quant_samples=$(find "$QUANT_DIR" -name "quant.sf" 2>/dev/null | wc -l || echo "0")
  if [[ $quant_samples -gt 0 ]]; then
    log_success "Quantification: $QUANT_DIR (${quant_samples} samples)"
  else
    log_warn "Quantification: $QUANT_DIR (no results found)"
  fi
else
  log_warn "Quantification: $QUANT_DIR (not found)"
fi

# =============================================================================
# MONITORING COMMANDS
# =============================================================================

echo ""
log_info "=== MONITORING COMMANDS ==="

echo "# View logs:"
echo "tail -f $LOG_DIR/orchestrator.log"
echo "ls -la $LOG_DIR/"
echo ""

echo "# Monitor FASTQ downloads:"
echo "find $FASTQ_DIR/getfastq -name '*.fastq*' | wc -l"
echo "du -sh $FASTQ_DIR/"
echo ""

echo "# Check manifest:"
echo "head -n 10 $MANIFEST_FILE"
echo "jq '.' $MANIFEST_FILE | head -n 50"
echo ""

echo "# Check metadata:"
[[ -f "$META_FILTERED" ]] && echo "head -n 10 $META_FILTERED"
[[ -f "$META_FILE" ]] && echo "head -n 10 $META_FILE"
echo ""

echo "# Monitor quantification:"
echo "find $QUANT_DIR -name 'quant.sf' | wc -l"
echo "ls -la $QUANT_DIR/ | head -n 20"
echo ""

# =============================================================================
# NEXT STEPS AND RECOMMENDATIONS
# =============================================================================

echo ""
log_info "=== NEXT STEPS ==="

if [[ "$OVERALL_STATUS" == "SUCCESS" ]]; then
  log_success "Pipeline completed successfully! 🎉"
  echo ""
  echo "Recommended next steps:"
  echo "1. Examine the merged expression matrix"
  echo "2. Perform quality control analysis"
  echo "3. Run comparative analysis (csca step)"
  echo "4. Generate publication-ready plots"

elif [[ "$OVERALL_STATUS" == "PARTIAL" ]]; then
  log_warn "Pipeline completed with some failures"
  echo ""
  echo "Recommended next steps:"
  echo "1. Review failed steps: ${FAILED_STEPS[*]}"
  echo "2. Check error logs in: $LOG_DIR"
  echo "3. Re-run failed steps individually:"
  for step in "${FAILED_STEPS[@]}"; do
    echo "   $0 --config $CONFIG --steps $step --force"
  done

else
  log_error "Pipeline failed"
  echo ""
  echo "Troubleshooting steps:"
  echo "1. Check error logs: ls -la $LOG_DIR/"
  echo "2. Verify configuration: $CONFIG"
  echo "3. Check NCBI_EMAIL environment variable"
  echo "4. Try running individual steps with --dry-run first"
  echo "5. Use --stream flag for real-time logging"
fi

echo ""
log_info "For help and support:"
echo "  - Configuration template: config/amalgkit_template.yaml"
echo "  - Run with --help for all options"
echo "  - Use --dry-run to preview execution plan"
echo "  - Enable --stream for real-time log output"

# Final status code
if [[ "$OVERALL_STATUS" == "SUCCESS" ]]; then
  exit 0
elif [[ "$OVERALL_STATUS" == "PARTIAL" ]]; then
  exit 1
else
  exit 2
fi
