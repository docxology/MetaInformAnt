#!/usr/bin/env bash
# Setup script for MetaInformAnt RNA-seq workflow on Linux
# This script installs all required dependencies

set -e  # Exit on error

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "================================"
echo "MetaInformAnt Setup - Linux"
echo "================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if running as root
if [ "$EUID" -eq 0 ]; then 
    echo -e "${RED}ERROR: Do not run this script as root${NC}"
    echo "Run as regular user. Script will prompt for sudo when needed."
    exit 1
fi

echo "Step 1: Checking Python and uv..."
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}ERROR: python3 not found${NC}"
    exit 1
fi
PYTHON_VERSION=$(python3 --version | awk '{print $2}')
echo -e "${GREEN}✓${NC} Python $PYTHON_VERSION"

if ! command -v uv &> /dev/null; then
    echo -e "${YELLOW}Installing uv...${NC}"
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.local/bin:$PATH"
fi
echo -e "${GREEN}✓${NC} uv $(uv --version | awk '{print $2}')"

echo ""
echo "Step 2: Creating Python virtual environment..."
if [ ! -d ".venv" ]; then
    uv venv
    echo -e "${GREEN}✓${NC} Virtual environment created"
else
    echo -e "${GREEN}✓${NC} Virtual environment exists"
fi

echo ""
echo "Step 3: Installing Python packages..."
source .venv/bin/activate

echo "  Installing metainformant..."
uv pip install -e . > /dev/null 2>&1
echo -e "${GREEN}✓${NC} metainformant installed"

echo "  Installing amalgkit..."
uv pip install git+https://github.com/kfuku52/amalgkit > /dev/null 2>&1
echo -e "${GREEN}✓${NC} amalgkit installed"

echo ""
echo "Step 4: Installing system dependencies..."

# Check if we can use sudo
if sudo -n true 2>/dev/null; then
    CAN_SUDO=true
else
    echo -e "${YELLOW}Note: May need sudo password for system package installation${NC}"
    CAN_SUDO=false
fi

# Install SRA Toolkit
if ! command -v fasterq-dump &> /dev/null; then
    echo "  Installing SRA Toolkit..."
    if command -v apt-get &> /dev/null; then
        sudo apt-get update -qq
        sudo apt-get install -y -qq sra-toolkit > /dev/null 2>&1
        echo -e "${GREEN}✓${NC} SRA Toolkit installed"
    else
        echo -e "${YELLOW}⚠${NC} Could not install SRA Toolkit automatically"
        echo "  Manual installation required: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit"
    fi
else
    echo -e "${GREEN}✓${NC} SRA Toolkit already installed"
fi

# Install kallisto
if ! command -v kallisto &> /dev/null; then
    echo "  Installing kallisto..."
    if command -v apt-get &> /dev/null; then
        sudo apt-get install -y -qq kallisto > /dev/null 2>&1
        echo -e "${GREEN}✓${NC} kallisto installed"
    else
        echo -e "${YELLOW}⚠${NC} Could not install kallisto automatically"
        echo "  Manual installation required: https://pachterlab.github.io/kallisto/download"
    fi
else
    echo -e "${GREEN}✓${NC} kallisto already installed"
fi

echo ""
echo "Step 5: Verifying installation..."

# Verify Python packages
python3 -c "import metainformant" 2>/dev/null && echo -e "${GREEN}✓${NC} metainformant module loads" || echo -e "${RED}✗${NC} metainformant module failed"
command -v amalgkit &>/dev/null && echo -e "${GREEN}✓${NC} amalgkit CLI available" || echo -e "${RED}✗${NC} amalgkit CLI missing"

# Verify system tools
command -v fasterq-dump &>/dev/null && echo -e "${GREEN}✓${NC} fasterq-dump available" || echo -e "${YELLOW}⚠${NC} fasterq-dump missing (optional but recommended)"
command -v prefetch &>/dev/null && echo -e "${GREEN}✓${NC} prefetch available" || echo -e "${YELLOW}⚠${NC} prefetch missing (optional)"
command -v kallisto &>/dev/null && echo -e "${GREEN}✓${NC} kallisto available" || echo -e "${YELLOW}⚠${NC} kallisto missing (required for quantification)"

echo ""
echo "Step 6: Environment configuration..."

# Check NCBI_EMAIL
if [ -z "$NCBI_EMAIL" ]; then
    echo -e "${YELLOW}⚠${NC} NCBI_EMAIL not set"
    echo "  Add to ~/.bashrc or ~/.zshrc:"
    echo "    export NCBI_EMAIL=\"your.email@example.com\""
else
    echo -e "${GREEN}✓${NC} NCBI_EMAIL=$NCBI_EMAIL"
fi

echo ""
echo "================================"
echo "Setup Complete!"
echo "================================"
echo ""
echo "To use the workflow:"
echo "  1. Activate virtual environment: source .venv/bin/activate"
echo "  2. Run workflow: python3 scripts/rna/run_multi_species.py"
echo ""
echo "For detailed documentation, see:"
echo "  docs/rna/SETUP.md"
echo ""

