#!/usr/bin/env bash
set -euo pipefail

# METAINFORMANT repo setup using uv and Python 3.11+
# - Creates/uses .venv
# - Installs project (and dev deps) via uv
# - Optionally installs amalgkit via uv
# - Optionally sets NCBI_EMAIL for this shell session and writes a helper file

usage() {
  cat <<EOF
Usage: bash scripts/setup_uv.sh [--with-amalgkit] [--ncbi-email EMAIL] [--skip-tests]

Options:
  --with-amalgkit       Install AMALGKIT from GitHub (optional)
  --ncbi-email EMAIL    Export NCBI_EMAIL for this session and write to output/setup/ncbi_email.txt
  --skip-tests          Do not run repository tests during setup

This script is idempotent and safe to re-run.
EOF
}

WITH_AMALGKIT=0
NCBI_EMAIL=""
DEFAULT_EMAIL="DanielAriFriedman@gmail.com"
SKIP_TESTS=0

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
    --skip-tests)
      SKIP_TESTS=1
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

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

if [[ ! -f pyproject.toml ]]; then
  echo "pyproject.toml not found; run from repo root" >&2
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

echo "[2/5] Creating/using virtual environment (.venv)"
UV_VENV_CLEAR=1 uv venv .venv

echo "[3/5] Installing project and dev dependencies"
uv pip install -e .[dev]

if [[ "$WITH_AMALGKIT" -eq 1 ]]; then
  echo "[4/5] Installing AMALGKIT via uv"
  uv pip install "git+https://github.com/kfuku52/amalgkit"
else
  echo "[4/5] Skipping AMALGKIT install (use --with-amalgkit to enable)"
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

echo "Verifying environment:"
uv run python -V
uv run metainformant --help >/dev/null 2>&1 && echo "metainformant CLI OK"

if [[ "$SKIP_TESTS" -eq 0 ]]; then
  echo "Optionally running tests (short):"
  uv run pytest -q || true
else
  echo "Skipping tests as requested (--skip-tests)"
fi

echo "[Extra] Ensuring MUSCLE CLI availability (for MSA tests)"
if ! command -v muscle >/dev/null 2>&1; then
  echo "MUSCLE not found in PATH. Installing lightweight fallback shim into .venv/bin/muscle"
  cat > .venv/bin/muscle <<'PYSHIM'
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
  chmod +x .venv/bin/muscle
  echo "Installed MUSCLE shim at .venv/bin/muscle"
else
  echo "MUSCLE found in PATH"
fi

echo "Setup complete. You can run commands with 'uv run <cmd>' or activate: source .venv/bin/activate"


