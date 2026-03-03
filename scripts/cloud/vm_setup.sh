#!/usr/bin/env bash
# VM native setup and pipeline launch
set -euo pipefail
LOGFILE="/opt/MetaInformAnt/vm_setup.log"
exec > >(tee -a "$LOGFILE") 2>&1
echo "=== VM Native Setup — $(date -Iseconds) ==="

# Install fastp
if ! command -v fastp &>/dev/null; then
    echo "▸ Installing fastp..."
    wget -q -O /usr/local/bin/fastp https://github.com/OpenGene/fastp/releases/download/v0.24.0/fastp
    chmod +x /usr/local/bin/fastp
    echo "  ✓ fastp $(fastp --version 2>&1)"
fi

# Install SRA Toolkit
if ! command -v fasterq-dump &>/dev/null; then
    echo "▸ Installing SRA Toolkit..."
    cd /tmp
    wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-ubuntu64.tar.gz
    tar xzf sratoolkit.3.2.1-ubuntu64.tar.gz
    cp sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump /usr/local/bin/
    cp sratoolkit.3.2.1-ubuntu64/bin/prefetch /usr/local/bin/
    cp sratoolkit.3.2.1-ubuntu64/bin/vdb-config /usr/local/bin/
    rm -rf sratoolkit*
    echo "  ✓ fasterq-dump installed"
fi

# Install Python deps
echo "▸ Installing Python dependencies..."
cd /opt/MetaInformAnt
pip3 install --break-system-packages -e "." 2>&1 | tail -3
pip3 install --break-system-packages amalgkit 2>&1 | tail -3
echo "  ✓ Python deps installed"

# Verify all tools
echo ""
echo "=== Tool Verification ==="
echo "kallisto: $(kallisto version 2>&1)"
echo "fastp: $(fastp --version 2>&1)"
echo "fasterq-dump: $(fasterq-dump --version 2>&1 | head -1)"
echo "python3: $(python3 --version)"
echo "amalgkit: $(pip3 show amalgkit 2>/dev/null | grep Version)"

# Create output directory
mkdir -p /opt/MetaInformAnt/output/amalgkit

# Start pipeline
echo ""
echo "=== Starting Pipeline — $(date -Iseconds) ==="
echo "Workers: 24 | Threads: 32 | Max GB: 20.0"
cd /opt/MetaInformAnt
nohup python3 scripts/rna/run_all_species.py --max-gb 20.0 --workers 24 --threads 32 > /opt/MetaInformAnt/output/amalgkit/pipeline.log 2>&1 &
echo "Pipeline PID: $!"
echo "✓ Pipeline started!"
echo "Monitor: tail -f /opt/MetaInformAnt/output/amalgkit/pipeline.log"
