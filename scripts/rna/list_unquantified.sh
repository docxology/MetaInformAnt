#!/bin/bash
# List unquantified samples that need processing
# This script identifies samples with downloaded SRA files but no quantification output

set -euo pipefail

# Determine paths relative to repository root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
AMALGKIT_DIR="${REPO_ROOT}/output/amalgkit"

# Color output
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m'

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_sample() {
    echo -e "${CYAN}  - $1${NC}"
}

# Function to list unquantified samples for a species
list_unquantified() {
    local species=$1
    local fastq_dir="${AMALGKIT_DIR}/${species}/fastq/getfastq"
    
    if [ ! -d "${fastq_dir}" ]; then
        return
    fi
    
    local samples=()
    local total_size=0
    
    for sample_dir in "${fastq_dir}"/*; do
        if [ ! -d "${sample_dir}" ]; then
            continue
        fi
        
        local sample=$(basename "${sample_dir}")
        local sra_file="${sample_dir}/${sample}.sra"
        local quant_file="${AMALGKIT_DIR}/${species}/quant/${sample}/abundance.tsv"
        
        if [ -f "${sra_file}" ] && [ ! -f "${quant_file}" ]; then
            local size=$(du -sh "${sample_dir}" | cut -f1)
            samples+=("${sample}")
            log_sample "${sample} (${size})"
            
            local size_kb=$(du -sk "${sample_dir}" | cut -f1)
            total_size=$((total_size + size_kb))
        fi
    done
    
    if [ ${#samples[@]} -gt 0 ]; then
        local size_gb=$(echo "scale=2; ${total_size}/1024/1024" | bc)
        echo -e "${YELLOW}  Total: ${#samples[@]} samples, ${size_gb}G${NC}"
        
        # Write sample list to output directory
        printf '%s\n' "${samples[@]}" > "${AMALGKIT_DIR}/${species}_unquantified.txt"
        echo -e "  Saved to: output/amalgkit/${species}_unquantified.txt"
    else
        echo "  (none)"
    fi
}

# Main execution
main() {
    log_info "Unquantified Samples Report"
    echo "======================================"
    echo
    
    for species in pbarbatus cfloridanus mpharaonis sinvicta; do
        if [ -d "${AMALGKIT_DIR}/${species}" ]; then
            echo "${species}:"
            list_unquantified "${species}"
            echo
        fi
    done
    
    log_info "Sample list files saved to: output/amalgkit/"
    echo
    log_info "Next steps:"
    echo "  1. Review output/amalgkit/*_unquantified.txt files for samples needing quantification"
    echo "  2. Run quantification workflow on these samples"
    echo "  3. After successful quantification, run scripts/rna/cleanup_quantified_sra.sh --execute"
}

main

