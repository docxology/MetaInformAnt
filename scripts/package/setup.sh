#!/usr/bin/env bash
set -euo pipefail

# METAINFORMANT repo setup using uv and Python 3.11+
# - Creates/uses .venv (or /tmp/metainformant_venv on FAT filesystems)
# - Installs project with dev and scientific dependencies via uv
# - Installs amalgkit via uv (default, use --skip-amalgkit to disable)
# - Optionally sets NCBI_EMAIL for this shell session and writes a helper file
# - Automatically detects FAT filesystems and configures UV cache accordingly
# - Supports --dev flag for development environment setup

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

usage() {
  cat <<EOF
Usage: bash scripts/setup.sh [--skip-amalgkit] [--with-deps] [--with-all] [--with-scraping] [--ncbi-email EMAIL] [--skip-tests] [--dev]

Options:
  --skip-amalgkit       Skip AMALGKIT installation (installed by default for RNA workflows)
  --with-deps           Install external CLI deps (seqkit, sra-tools, kallisto, R); uv pip install parallel-fastq-dump
  --with-all            Install all optional dependencies (database, networks, scraping, etc.) - scientific deps installed by default
  --with-scraping       Install scraping dependencies (cloudscraper) for web scraping functionality
  --ncbi-email EMAIL    Export NCBI_EMAIL for this session and write to output/setup/ncbi_email.txt
  --skip-tests          Do not run repository tests during setup
  --dev                 Setup development environment (pre-commit, docs, output dirs, wrapper scripts)

Default behavior: Installs dev, scientific dependencies (scipy, scikit-learn, etc.), and amalgkit

This script is idempotent and safe to re-run.
EOF
}

WITH_AMALGKIT=1  # Install amalgkit by default
NCBI_EMAIL=""
DEFAULT_EMAIL="DanielAriFriedman@gmail.com"
SKIP_TESTS=0
WITH_DEPS=0
WITH_SCIENTIFIC=0
WITH_ALL=0
WITH_SCRAPING=0
WITH_DEV=0  # Development environment setup

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-amalgkit)
      WITH_AMALGKIT=0
      shift
      ;;
    --with-amalgkit)
      # Keep for backward compatibility, but it's now default
      WITH_AMALGKIT=1
      shift
      ;;
    --with-deps)
      WITH_DEPS=1
      shift
      ;;
    --with-scientific)
      # Keep for backward compatibility, but it's now default
      WITH_SCIENTIFIC=1
      shift
      ;;
    --with-all)
      WITH_ALL=1
      shift
      ;;
    --with-scraping)
      WITH_SCRAPING=1
      shift
      ;;
    --ncbi-email)
      NCBI_EMAIL="${2:-}"
      shift 2 || true
      ;;
    --skip-tests)
      SKIP_TESTS=1
      shift
      ;;
    --dev)
      WITH_DEV=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 2
      ;;
  esac
done

# Get the absolute path of the script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Go up from scripts/package to repo root
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$ROOT_DIR"

if [[ ! -f pyproject.toml ]]; then
  echo "pyproject.toml not found; run from repo root" >&2
  echo "Current directory: $(pwd)" >&2
  exit 2
fi

echo "[1/5] Ensuring uv is available..."
if ! command -v uv >/dev/null 2>&1; then
  cat >&2 <<'EOM'
uv not found on PATH.
Install uv using one of the following and re-run this script:
  - macOS (Homebrew):   brew install uv
  - Official installer: curl -LsSf https://astral.sh/uv/install.sh | sh
See: https://github.com/astral-sh/uv
EOM
  exit 3
fi

echo "[1.5/5] Detecting filesystem type and configuring UV cache..."
# Detect filesystem type and set UV_CACHE_DIR if needed
# Use Python to detect filesystem (works even if package not installed yet)
# Temporarily disable exit on error to capture Python script exit code
set +e
FS_DETECT=$(python3 <<'PYTHON'
import subprocess
import sys
from pathlib import Path

def detect_fs_type(path):
    """Detect filesystem type using df -T."""
    try:
        result = subprocess.run(
            ["df", "-T", str(path)],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            lines = result.stdout.strip().split("\n")
            if len(lines) > 1:
                parts = lines[1].split()
                if len(parts) >= 2:
                    return parts[1].lower()
    except Exception:
        pass
    return "unknown"

def supports_symlinks(fs_type):
    """Check if filesystem supports symlinks."""
    no_symlink_fs = {"exfat", "fat32", "fat", "vfat", "msdos"}
    return fs_type not in no_symlink_fs

repo_root = Path.cwd()
fs_type = detect_fs_type(repo_root)
has_symlinks = supports_symlinks(fs_type)

print(f"Filesystem type: {fs_type}")
print(f"Symlink support: {has_symlinks}")

if not has_symlinks:
    cache_dir = "/tmp/uv-cache"
    venv_dir = "/tmp/metainformant_venv"
    print(f"FAT filesystem detected - using {cache_dir} for UV cache")
    print(f"Virtual environment will be created at {venv_dir}")
    sys.exit(1)  # Signal to use alternative locations
else:
    cache_dir = ".uv-cache"
    venv_dir = ".venv"
    print(f"Using standard locations: cache={cache_dir}, venv={venv_dir}")
    sys.exit(0)
PYTHON
)

FS_EXIT_CODE=$?
set -e

if [[ $FS_EXIT_CODE -eq 1 ]]; then
  # FAT filesystem detected - use /tmp locations
  export UV_CACHE_DIR="/tmp/uv-cache"
  VENV_DIR="/tmp/metainformant_venv"
  echo "  → FAT filesystem detected - using /tmp/uv-cache for UV cache"
  echo "  → Virtual environment will be created at $VENV_DIR"
else
  # Standard filesystem - use repo locations
  VENV_DIR=".venv"
  echo "  → Standard filesystem - using .uv-cache for UV cache"
fi

# Create cache directory if needed
if [[ -n "$UV_CACHE_DIR" ]]; then
  mkdir -p "$UV_CACHE_DIR"
  echo "  → UV cache directory: $UV_CACHE_DIR"
fi

echo "[2/5] Creating/using virtual environment ($VENV_DIR)"
if [[ "$VENV_DIR" == "/tmp/metainformant_venv" ]]; then
  # For FAT filesystem, try to create venv in /tmp
  if [[ -d "$VENV_DIR" ]]; then
    echo "  → Using existing venv at $VENV_DIR"
  else
    UV_VENV_CLEAR=1 uv venv "$VENV_DIR" || {
      echo "Warning: Failed to create venv at $VENV_DIR, trying .venv..." >&2
      UV_VENV_CLEAR=1 uv venv .venv
      VENV_DIR=".venv"
    }
  fi
else
  UV_VENV_CLEAR=1 uv venv "$VENV_DIR"
fi

echo "[3/5] Installing project and dependencies"
# Build dependency list - install scientific deps by default
DEPS="dev,scientific"
if [[ "$WITH_ALL" -eq 1 ]]; then
  DEPS="${DEPS},all"
elif [[ "$WITH_SCRAPING" -eq 1 ]]; then
  DEPS="${DEPS},scraping"
fi

# Use venv Python explicitly if using alternative location
if [[ "$VENV_DIR" == "/tmp/metainformant_venv" ]]; then
  uv pip install -e ".[${DEPS}]" --python "$VENV_DIR/bin/python3"
else
  uv pip install -e ".[${DEPS}]"
fi

if [[ "$WITH_ALL" -eq 1 ]]; then
  echo "  → Installed all optional dependencies (scientific, database, networks, scraping, etc.)"
elif [[ "$WITH_SCRAPING" -eq 1 ]]; then
  echo "  → Installed scientific dependencies (scipy, scikit-learn, etc.)"
  echo "  → Installed scraping dependencies (cloudscraper)"
  echo "  → Use --with-all to install additional optional dependencies (database, networks, etc.)"
else
  echo "  → Installed scientific dependencies (scipy, scikit-learn, etc.)"
  echo "  → Use --with-scraping to install scraping dependencies (cloudscraper)"
  echo "  → Use --with-all to install all optional dependencies (database, networks, scraping, etc.)"
fi

if [[ "$WITH_AMALGKIT" -eq 1 ]]; then
  echo "[4/5] Installing AMALGKIT via uv"
  if [[ "$VENV_DIR" == "/tmp/metainformant_venv" ]]; then
    uv pip install "git+https://github.com/kfuku52/amalgkit" --python "$VENV_DIR/bin/python3"
  else
    uv pip install "git+https://github.com/kfuku52/amalgkit"
  fi
  # Verify amalgkit installation
  if [[ "$VENV_DIR" == "/tmp/metainformant_venv" ]]; then
    if "$VENV_DIR/bin/python3" -c "import amalgkit" 2>/dev/null; then
      echo "  → ✓ AMALGKIT installed successfully"
    else
      echo "  → ⚠ AMALGKIT installation may have failed - check above for errors"
    fi
  else
    if uv run python -c "import amalgkit" 2>/dev/null; then
      echo "  → ✓ AMALGKIT installed successfully"
    else
      echo "  → ⚠ AMALGKIT installation may have failed - check above for errors"
    fi
  fi
else
  echo "[4/5] Skipping AMALGKIT install (use --skip-amalgkit to disable, or omit flag to install)"
fi

mkdir -p output/setup

# Decide which email to use
if [[ -z "$NCBI_EMAIL" ]]; then
  if [[ -f output/setup/ncbi_email.txt ]]; then
    NCBI_EMAIL="$(cat output/setup/ncbi_email.txt)"
  else
    NCBI_EMAIL="$DEFAULT_EMAIL"
    echo "[5/5] Using default NCBI_EMAIL: $NCBI_EMAIL"
  fi
fi

if [[ -n "$NCBI_EMAIL" ]]; then
  echo "[5/5] Exporting NCBI_EMAIL for this session and recording to output/setup/ncbi_email.txt"
  export NCBI_EMAIL="$NCBI_EMAIL"
  printf "%s\n" "$NCBI_EMAIL" > output/setup/ncbi_email.txt
else
  echo "[5/5] NCBI_EMAIL not provided; you can pass --ncbi-email to set it"
fi

# Optional: install external CLI dependencies
if [[ "$WITH_DEPS" -eq 1 ]]; then
  echo "Installing external CLI dependencies (requires Homebrew on macOS)..."
  if command -v brew >/dev/null 2>&1; then
    # Ensure useful taps are present for bio tools
    brew tap brewsci/bio || true
    # Install individually to avoid one missing formula aborting the rest
    for pkg in wget curl pigz parallel seqkit samtools kallisto r; do
      echo "brew install $pkg (best-effort)"
      brew install "$pkg" || true
    done
    # Salmon may live under brewsci/bio
    echo "Attempting to install salmon (best-effort)"
    brew install salmon || brew install brewsci/bio/salmon || true
    # SRA Toolkit may be packaged as sratoolkit; try formula and cask
    echo "Attempting to install sratoolkit (best-effort)"
    brew install sratoolkit || brew install --cask sratoolkit || true
    # Ensure Homebrew bin dir is on PATH for this shell
    if [[ -d "/opt/homebrew/bin" ]]; then
      export PATH="/opt/homebrew/bin:$PATH"
    fi
  else
    echo "Homebrew not found; please install wget, curl, pigz, parallel, seqkit, samtools, kallisto, salmon, sratoolkit, and R manually." >&2
  fi
  # Install parallel-fastq-dump as a Python console script into the venv
  if [[ "$VENV_DIR" == "/tmp/metainformant_venv" ]]; then
    uv pip install parallel-fastq-dump --python "$VENV_DIR/bin/python3" || true
  else
    uv pip install parallel-fastq-dump || true
  fi
fi

echo "Verifying environment:"
if [[ "$VENV_DIR" == "/tmp/metainformant_venv" ]]; then
  "$VENV_DIR/bin/python3" -V
  "$VENV_DIR/bin/python3" -m metainformant --help >/dev/null 2>&1 && echo "metainformant CLI OK"
else
  uv run python -V
  uv run metainformant --help >/dev/null 2>&1 && echo "metainformant CLI OK"
fi

if [[ "$SKIP_TESTS" -eq 0 ]]; then
  echo "Running tests:"
  echo "  → This may take a few minutes..."
  echo "  → Tests are running - you'll see test names as they execute..."
  echo "  → Progress shown in real-time below..."
  echo ""
  if [[ "$VENV_DIR" == "/tmp/metainformant_venv" ]]; then
    PYTEST_CMD="$VENV_DIR/bin/pytest"
  else
    PYTEST_CMD="uv run pytest"
  fi
  
  # Count total tests first for progress indication
  echo "  → Collecting test files..."
  TEST_COUNT=$($PYTEST_CMD --collect-only -q 2>/dev/null | grep -c "test session starts\|<" || echo "0")
  if [[ "$TEST_COUNT" -gt 0 ]]; then
    echo "  → Found $TEST_COUNT test(s) to run"
  fi
  echo ""
  
  # Run tests with real-time output showing test progress
  # Use -v for verbose (shows test names), --tb=short for cleaner errors
  # Show output in real-time while filtering import errors
  if command -v stdbuf >/dev/null 2>&1; then
    # Use stdbuf for unbuffered output if available
    $PYTEST_CMD -v --tb=short --durations=0 2>&1 | \
      stdbuf -oL -eL grep -v --line-buffered \
        -e "ModuleNotFoundError: No module named 'scipy'" \
        -e "ModuleNotFoundError: No module named 'psycopg2'" \
        -e "^ERROR collecting tests/test_core_db.py" \
        -e "^ERROR collecting tests/test_gwas" | \
      tee /tmp/pytest_$$.log || true
  else
    # Fallback: use unbuffered Python if stdbuf not available
    $PYTEST_CMD -v --tb=short --durations=0 2>&1 | \
      python3 -u -c "import sys; [sys.stdout.write(l) for l in sys.stdin if not any(x in l for x in ['ModuleNotFoundError: No module named \\'scipy\\'', 'ModuleNotFoundError: No module named \\'psycopg2\\'', 'ERROR collecting tests/test_core_db.py', 'ERROR collecting tests/test_gwas'])]" | \
      tee /tmp/pytest_$$.log || true
  fi
  
  TEST_OUTPUT=$(cat /tmp/pytest_$$.log 2>/dev/null || echo "")
  rm -f /tmp/pytest_$$.log
  
  # Show test results summary
  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "Test Results Summary:"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  
  # Extract summary lines (pytest prints these at the end)
  if [[ -n "$TEST_OUTPUT" ]]; then
    # Look for pytest summary lines (multiple formats)
    # Format 1: "X passed in Y.YYs"
    # Format 2: "X passed, Y failed"
    # Format 3: "ERRORS" section
    SUMMARY=$(echo "$TEST_OUTPUT" | grep -E "([0-9]+ (passed|failed|error|warnings|skipped)|passed in|failed in|error in)" | tail -5)
    if [[ -n "$SUMMARY" ]]; then
      echo "$SUMMARY"
    else
      # Fallback: show last few lines if no summary found
      echo "$TEST_OUTPUT" | tail -10
    fi
  else
    echo "  → Test output shown above (real-time progress)"
  fi
  
  # Check if any optional deps are missing and show helpful message
  if [[ "$WITH_ALL" -eq 0 ]] && [[ -n "$TEST_OUTPUT" ]]; then
    if echo "$TEST_OUTPUT" | grep -q "ModuleNotFoundError: No module named 'psycopg2'"; then
      echo ""
      echo "  ℹ️  Some tests require database dependencies (psycopg2, etc.)"
      echo "     Install with: bash scripts/package/setup_uv.sh --with-all"
    fi
  fi
  
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
else
  echo "Skipping tests as requested (--skip-tests)"
fi

echo "[Extra] Ensuring MUSCLE CLI availability (for MSA tests)"
MUSCLE_PATH="$VENV_DIR/bin/muscle"
if ! command -v muscle >/dev/null 2>&1; then
  echo "MUSCLE not found in PATH. Installing lightweight fallback shim into $MUSCLE_PATH"
  cat > "$MUSCLE_PATH" <<'PYSHIM'
#!/usr/bin/env python3
import sys, argparse
from pathlib import Path
from Bio import SeqIO
from metainformant.dna.msa import align_msa

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-align', dest='align', help='input FASTA')
    p.add_argument('-output', dest='output', help='output FASTA')
    # Accept common no-op flags to be tolerant
    p.add_argument('-quiet', action='store_true')
    p.add_argument('-threads', type=int, default=1)
    args, _ = p.parse_known_args()
    if not args.align or not args.output:
        sys.stderr.write('muscle shim error: require -align <in> -output <out>\n')
        sys.exit(2)
    id_to_seq = {}
    for rec in SeqIO.parse(args.align, 'fasta'):
        id_to_seq[rec.id] = str(rec.seq)
    aligned = align_msa(id_to_seq)
    with open(args.output, 'w') as fh:
        for k, v in aligned.items():
            fh.write(f'>{k}\n{v}\n')
    return 0

if __name__ == '__main__':
    sys.exit(main())
PYSHIM
  chmod +x "$MUSCLE_PATH"
  echo "Installed MUSCLE shim at $MUSCLE_PATH"
else
  echo "MUSCLE found in PATH"
fi

# Development environment setup (if requested)
if [[ "$WITH_DEV" -eq 1 ]]; then
  echo ""
  echo "[DEV] Setting up development environment..."

  # Install pre-commit hooks
  echo "  Installing pre-commit hooks..."
  uv run pre-commit install

  # Create output directory structure following .cursorrules
  echo "  Creating output directory structure..."
  mkdir -p output/{coverage_html,benchmarks,profiles,docs,test_results}
  mkdir -p output/{dna,rna,protein,math,networks,singlecell,visualization}
  mkdir -p output/{simulation,ontology,quality,multiomics,epigenome}

  # Set up documentation build environment
  echo "  Setting up documentation environment..."
  mkdir -p docs/_build docs/_static docs/_templates

  echo "  ✅ Development environment setup complete"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ Setup complete!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if [[ "$VENV_DIR" == "/tmp/metainformant_venv" ]]; then
  echo "Virtual environment: $VENV_DIR"
  echo "UV cache directory: $UV_CACHE_DIR"
  echo ""
  echo "Usage:"
  echo "  $VENV_DIR/bin/python3 -m metainformant --help"
  echo "  source $VENV_DIR/bin/activate"
else
  echo "Virtual environment: $VENV_DIR"
  echo ""
  echo "Usage:"
  echo "  uv run metainformant --help"
  echo "  source $VENV_DIR/bin/activate"
fi

if [[ "$WITH_ALL" -eq 0 ]]; then
  echo ""
  echo "💡 Tip: Install optional dependencies for extended functionality:"
  if [[ "$WITH_SCRAPING" -eq 0 ]]; then
    echo "  - Scraping (cloudscraper): bash scripts/package/setup.sh --with-scraping"
  fi
  echo "  - All dependencies: bash scripts/package/setup.sh --with-all"
  if [[ "$WITH_DEV" -eq 0 ]]; then
    echo "  - Development environment: bash scripts/package/setup.sh --dev"
  fi
fi
