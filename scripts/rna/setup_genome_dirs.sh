#!/usr/bin/env bash
# setup_genome_dirs.sh – Create output directory structure for large-scale amalgkit RNA-seq runs
#
# Usage:
#   bash scripts/rna/setup_genome_dirs.sh                        # create all known species dirs
#   bash scripts/rna/setup_genome_dirs.sh --species "apis_mellifera pogonomyrmex_barbatus"
#
# Creates incremental storage layout so quant data accumulates safely across runs:
#   output/amalgkit/{species}/work/{metadata,index,quant}/   ← PERMANENT (never delete)
#   output/amalgkit/{species}/fastq/                          ← TEMPORARY (auto-deleted after quant)
#   output/amalgkit/{species}/{merged,curate,cstmm,csca,logs}/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT_DIR"

GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

info() { echo -e "${BLUE}ℹ️${NC}  $*"; }
ok()   { echo -e "${GREEN}✅${NC} $*"; }

# Default: all species with existing configs
CUSTOM_SPECIES=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --species)
      CUSTOM_SPECIES="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: bash $0 [--species \"species1 species2\"]"
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

# Build species list from config filenames if not specified
if [[ -n "$CUSTOM_SPECIES" ]]; then
  SPECIES_LIST=($CUSTOM_SPECIES)
else
  # Derive species names from config files (e.g., amalgkit_apis_mellifera.yaml → apis_mellifera)
  SPECIES_LIST=()
  while IFS= read -r cfg; do
    base=$(basename "$cfg" .yaml)
    # Strip amalgkit_ prefix, skip template/test/cross
    species="${base#amalgkit_}"
    case "$species" in
      template|test|cross_species|faq) continue ;;
      *) SPECIES_LIST+=("$species") ;;
    esac
  done < <(find config/amalgkit -name "amalgkit_*.yaml" | sort)
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " METAINFORMANT – RNA-seq Output Directory Setup"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
info "Creating directories for ${#SPECIES_LIST[@]} species under output/amalgkit/"
echo ""

for species in "${SPECIES_LIST[@]}"; do
  echo "  → $species"

  # Permanent directories: quant data accumulates here incrementally
  PERM_DIRS=(
    "output/amalgkit/$species/work/metadata"    # NCBI SRA metadata
    "output/amalgkit/$species/work/index"       # kallisto index (built once)
    "output/amalgkit/$species/work/quant"       # per-sample abundance.tsv (~2 MB each)
    "output/amalgkit/$species/work/logs"        # workflow logs
  )

  # Temporary directories: FASTQs are downloaded here then deleted after quantification
  TEMP_DIRS=(
    "output/amalgkit/$species/fastq"            # temp FASTQs (auto-deleted)
  )

  # Post-processing output directories
  POST_DIRS=(
    "output/amalgkit/$species/merged"           # merged expression matrix
    "output/amalgkit/$species/cstmm"            # cross-species TMM normalization
    "output/amalgkit/$species/curate"           # outlier removal + QC plots
    "output/amalgkit/$species/csca"             # cross-species correlation
    "output/amalgkit/$species/logs"             # species-level logs
  )

  for d in "${PERM_DIRS[@]}" "${TEMP_DIRS[@]}" "${POST_DIRS[@]}"; do
    mkdir -p "$d"
  done
done

echo ""
ok "Directory structure created for ${#SPECIES_LIST[@]} species"

# Create a .gitignore to prevent accidental commit of large data files
cat > output/amalgkit/.gitignore << 'EOF'
# Large data directories – do NOT commit to git
*/fastq/
*/work/quant/
*/work/index/*.idx
*/merged/*.tsv.gz
*/cstmm/
*/curate/*/plots/
# Keep metadata (small, useful to track)
!*/work/metadata/
EOF

ok "Created output/amalgkit/.gitignore to protect large data files"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " Storage Layout Reminder"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  PERMANENT (never delete – incremental quant data):"
echo "    output/amalgkit/*/work/quant/    ← abundance.tsv per sample"
echo "    output/amalgkit/*/work/index/    ← kallisto index"
echo "    output/amalgkit/*/work/metadata/ ← SRA metadata"
echo ""
echo "  TEMPORARY (auto-deleted after quantification):"
echo "    output/amalgkit/*/fastq/         ← downloaded FASTQs"
echo ""
echo "  SAFE TO DELETE (regenerable):"
echo "    output/amalgkit/*/merged/        ← re-merge from quant/"
echo "    output/amalgkit/*/curate/        ← re-run curate step"
echo ""
