#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 WHEEL PYTHON_VERSION" >&2
  exit 2
fi

wheel_path="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
python_version="$2"
tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/metainformant-wheel-smoke.XXXXXX")"
trap 'rm -rf "$tmp_dir"' EXIT

uv venv --python "$python_version" "$tmp_dir/venv"
uv pip install --python "$tmp_dir/venv/bin/python" "$wheel_path"

"$tmp_dir/venv/bin/python" - <<'PY'
import metainformant
from importlib.util import find_spec
from metainformant.core.io import paths
from metainformant.core.utils import config
from metainformant.dna.sequence import core

assert metainformant.__version__ == "0.4.0"
assert callable(paths.get_project_root)
assert callable(config.load_mapping_from_file)
assert callable(core.gc_content)
for removed in ("metainformant.core.config", "metainformant.core.paths", "metainformant.dna.sequences"):
    assert find_spec(removed) is None, f"removed compatibility module still installed: {removed}"
print(f"wheel smoke import OK: {metainformant.__version__}")
PY

"$tmp_dir/venv/bin/python" -m metainformant --help >/dev/null
echo "wheel smoke passed for Python ${python_version}: ${wheel_path}"
