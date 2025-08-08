#!/usr/bin/env bash
set -euo pipefail

# METAINFORMANT repo setup using uv and Python 3.11+
# - Creates/uses .venv
# - Installs project (and dev deps)
# - Optionally installs amalgkit
# - Optionally sets NCBI_EMAIL for this shell session and writes a helper file

usage() {
  cat <<EOF
Usage: bash scripts/setup_uv.sh [--with-amalgkit] [--ncbi-email EMAIL]

Options:
  --with-amalgkit       Install AMALGKIT from GitHub (optional)
  --ncbi-email EMAIL    Export NCBI_EMAIL for this session and write to output/setup/ncbi_email.txt

This script is idempotent and safe to re-run.
EOF
}

WITH_AMALGKIT=0
NCBI_EMAIL=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-amalgkit)
      WITH_AMALGKIT=1
      shift
      ;;
    --ncbi-email)
      NCBI_EMAIL="${2:-}"
      shift 2 || true
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

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

if [[ ! -f pyproject.toml ]]; then
  echo "pyproject.toml not found; run from repo root" >&2
  exit 2
fi

echo "[1/5] Ensuring uv is available..."
if ! command -v uv >/dev/null 2>&1; then
  echo "uv not found on PATH. Attempting pip install --user uv"
  python3 -m pip install --user uv >/dev/null 2>&1 || true
  if ! command -v uv >/dev/null 2>&1; then
    echo "uv still not found. Please install uv (https://github.com/astral-sh/uv) and re-run." >&2
    exit 3
  fi
fi

echo "[2/5] Creating/using virtual environment (.venv)"
uv venv .venv
source .venv/bin/activate

echo "[3/5] Installing project and dev dependencies"
uv pip install -e .[dev]

if [[ "$WITH_AMALGKIT" -eq 1 ]]; then
  echo "[4/5] Installing AMALGKIT"
  pip install --no-input --no-warn-script-location "git+https://github.com/kfuku52/amalgkit"
else
  echo "[4/5] Skipping AMALGKIT install (use --with-amalgkit to enable)"
fi

mkdir -p output/setup

if [[ -n "$NCBI_EMAIL" ]]; then
  echo "[5/5] Exporting NCBI_EMAIL for this session and recording to output/setup/ncbi_email.txt"
  export NCBI_EMAIL="$NCBI_EMAIL"
  printf "%s\n" "$NCBI_EMAIL" > output/setup/ncbi_email.txt
else
  echo "[5/5] NCBI_EMAIL not provided; you can pass --ncbi-email to set it"
fi

echo "Verifying environment:"
python -V
metainformant --help >/dev/null 2>&1 && echo "metainformant CLI OK"

echo "Optionally running tests (short):"
pytest -q || true

echo "Setup complete. To activate later: source .venv/bin/activate"


