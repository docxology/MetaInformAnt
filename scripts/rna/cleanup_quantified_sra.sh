#!/bin/bash
# Cleanup script for quantified SRA files in AmalgKit output
# This script safely removes SRA files only after verifying quantification is complete

set -euo pipefail

# Determine paths relative to repository root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
AMALGKIT_DIR="${REPO_ROOT}/output/amalgkit"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if sample is quantified
is_quantified() {
    local species=$1
    local sample=$2
    local quant_dir="${AMALGKIT_DIR}/${species}/quant/${sample}"
    
    if [ -f "${quant_dir}/abundance.tsv" ] && [ -f "${quant_dir}/run_info.json" ]; then
        return 0
    else
        return 1
    fi
}

# Function to cleanup SRA files for a species
cleanup_species() {
    local species=$1
    local fastq_dir="${AMALGKIT_DIR}/${species}/fastq/getfastq"
    
    if [ ! -d "${fastq_dir}" ]; then
        log_warn "No fastq directory found for ${species}"
        return
    fi
    
    log_info "Processing ${species}..."
    
    local deleted_count=0
    local deleted_size=0
    local skipped_count=0
    
    for sample_dir in "${fastq_dir}"/*; do
        if [ ! -d "${sample_dir}" ]; then
            continue
        fi
        
        local sample=$(basename "${sample_dir}")
        local sra_file="${sample_dir}/${sample}.sra"
        
        if [ ! -f "${sra_file}" ]; then
            continue
        fi
        
        if is_quantified "${species}" "${sample}"; then
            local size=$(du -sk "${sra_file}" | cut -f1)
            log_info "  Deleting ${sample}.sra ($(echo "scale=2; ${size}/1024" | bc)M)"
            rm -f "${sra_file}"
            deleted_count=$((deleted_count + 1))
            deleted_size=$((deleted_size + size))
        else
            log_warn "  Skipping ${sample}.sra (not quantified)"
            skipped_count=$((skipped_count + 1))
        fi
    done
    
    if [ ${deleted_count} -gt 0 ]; then
        local size_gb=$(echo "scale=2; ${deleted_size}/1024/1024" | bc)
        log_info "  Deleted ${deleted_count} SRA files, reclaimed ${size_gb}G"
    fi
    
    if [ ${skipped_count} -gt 0 ]; then
        log_warn "  Skipped ${skipped_count} unquantified samples"
    fi
}

# Main execution
main() {
    log_info "Starting SRA cleanup for quantified samples"
    log_info "Working directory: ${AMALGKIT_DIR}"
    echo
    
    # Process each species
    for species in pbarbatus cfloridanus mpharaonis sinvicta; do
        if [ -d "${AMALGKIT_DIR}/${species}" ]; then
            cleanup_species "${species}"
            echo
        fi
    done
    
    log_info "Cleanup complete!"
    log_info "Generating updated storage report..."
    
    # Show updated sizes
    echo
    log_info "Updated storage usage:"
    cd "${AMALGKIT_DIR}"
    du -sh */ 2>/dev/null | sort -h
}

# Run with confirmation
if [ "$#" -eq 1 ] && [ "$1" = "--execute" ]; then
    main
else
    echo "This script will delete quantified SRA files to reclaim disk space."
    echo
    echo "P. barbatus: ~19.83 GB can be reclaimed from 18 quantified samples"
    echo "C. floridanus: No quantified samples yet (needs processing first)"
    echo
    echo "Run with --execute flag to proceed:"
    echo "  $0 --execute"
    exit 0
fi

