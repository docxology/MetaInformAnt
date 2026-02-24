#!/usr/bin/env bash
# run_all_species.sh – Run all amalgkit species workflows sequentially
#
# Usage:
#   bash scripts/rna/run_all_species.sh              # run all species
#   bash scripts/rna/run_all_species.sh --dry-run    # show what would run
#   bash scripts/rna/run_all_species.sh --status     # check status of all species
#
# Each species runs to completion before the next starts.
# Safe to interrupt and re-run — each workflow resumes from where it left off.
# Log for each species: output/amalgkit/logs/{species}_{timestamp}.log

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT_DIR"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

info()  { echo -e "${BLUE}[$(date '+%H:%M:%S')]${NC} $*"; }
ok()    { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✅${NC} $*"; }
warn()  { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠️${NC}  $*"; }
fail()  { echo -e "${RED}[$(date '+%H:%M:%S')] ❌${NC} $*"; }

DRY_RUN=0
STATUS_ONLY=0
SPECIES_FILTER=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)   DRY_RUN=1 ;;
    --status)    STATUS_ONLY=1 ;;
    --species)   SPECIES_FILTER="$2"; shift ;;
    -h|--help)
      echo "Usage: bash $0 [--dry-run] [--status] [--species <name>]"
      exit 0
      ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
  shift
done

# ── Collect configs (sorted, skip template/test/cross) ──────────────────────
CONFIGS=()
while IFS= read -r f; do
  base=$(basename "$f" .yaml)
  species="${base#amalgkit_}"
  # If filter specified, only run matching species
  if [[ -n "$SPECIES_FILTER" && "$species" != *"$SPECIES_FILTER"* ]]; then
    continue
  fi
  CONFIGS+=("$f")
done < <(ls config/amalgkit/amalgkit_*.yaml | sort | grep -v -E "template|_test\.yaml|cross_species")

TOTAL=${#CONFIGS[@]}

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " METAINFORMANT – Sequential Species Processor"
echo " Species to process: $TOTAL"
echo " Threads per species: 16 (chunk-size 16, per-sample concurrency)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Activate venv
if [[ -f .venv/bin/activate ]]; then
  source .venv/bin/activate
elif [[ -f /tmp/metainformant_venv/bin/activate ]]; then
  source /tmp/metainformant_venv/bin/activate
fi

# Status mode: just check everything
if [[ "$STATUS_ONLY" -eq 1 ]]; then
  for cfg in "${CONFIGS[@]}"; do
    species=$(basename "$cfg" .yaml | sed 's/amalgkit_//')
    echo ""
    info "Status: $species"
    python3 scripts/rna/run_workflow.py --config "$cfg" --status 2>&1 | tail -5 || true
  done
  exit 0
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo ""
  info "DRY RUN – would process:"; for cfg in "${CONFIGS[@]}"; do echo "  $(basename "$cfg")"; done
  exit 0
fi

# ── Main sequential loop ─────────────────────────────────────────────────────
PASSED=0
FAILED=0
FAILED_SPECIES=()
OVERALL_LOG="output/amalgkit/run_all_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p output/amalgkit
echo "Overall log: $OVERALL_LOG"
echo "Started: $(date)" | tee "$OVERALL_LOG"

for idx in "${!CONFIGS[@]}"; do
  cfg="${CONFIGS[$idx]}"
  species=$(basename "$cfg" .yaml | sed 's/amalgkit_//')
  num=$((idx + 1))

  echo "" | tee -a "$OVERALL_LOG"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$OVERALL_LOG"
  info "[$num/$TOTAL] Starting: $species" | tee -a "$OVERALL_LOG"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$OVERALL_LOG"

  # Per-species log file
  mkdir -p "output/amalgkit/$species/logs"
  SPECIES_LOG="output/amalgkit/$species/logs/workflow_$(date '+%Y%m%d_%H%M%S').log"

  START_TIME=$(date +%s)

  if python3 scripts/rna/run_workflow.py --config "$cfg" --stream --chunk-size 16 2>&1 | tee "$SPECIES_LOG" | tee -a "$OVERALL_LOG"; then
    END_TIME=$(date +%s)
    ELAPSED=$(( (END_TIME - START_TIME) / 60 ))
    ok "[$num/$TOTAL] Completed: $species (${ELAPSED}m)" | tee -a "$OVERALL_LOG"
    PASSED=$((PASSED + 1))
  else
    END_TIME=$(date +%s)
    ELAPSED=$(( (END_TIME - START_TIME) / 60 ))
    fail "[$num/$TOTAL] FAILED: $species (${ELAPSED}m) – see $SPECIES_LOG" | tee -a "$OVERALL_LOG"
    FAILED=$((FAILED + 1))
    FAILED_SPECIES+=("$species")
    warn "Continuing to next species..." | tee -a "$OVERALL_LOG"
  fi
done

# ── Summary ──────────────────────────────────────────────────────────────────
echo "" | tee -a "$OVERALL_LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$OVERALL_LOG"
echo " FINAL SUMMARY" | tee -a "$OVERALL_LOG"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" | tee -a "$OVERALL_LOG"
ok "Completed: $PASSED / $TOTAL species" | tee -a "$OVERALL_LOG"
if [[ "$FAILED" -gt 0 ]]; then
  fail "Failed: $FAILED species: ${FAILED_SPECIES[*]}" | tee -a "$OVERALL_LOG"
fi
echo "Finished: $(date)" | tee -a "$OVERALL_LOG"
echo "Overall log: $OVERALL_LOG" | tee -a "$OVERALL_LOG"

[[ "$FAILED" -eq 0 ]]
