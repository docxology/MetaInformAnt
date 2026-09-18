#!/bin/bash
# Blocking security and dependency checks for the MetaInformAnt release surface.
#
# Mirrors the blocking security-checks job in .github/workflows/security.yml so
# developers can reproduce the CI gate locally before pushing. Three gates:
#
#   1. bandit static analysis over src/metainformant (tests excluded)
#   2. pip-audit vulnerability scan over the uv.lock-derived requirements
#   3. grep audit for shell-execution and path-traversal patterns in
#      maintained src/ paths
#
# Gate policy: medium-and-higher severity bandit findings, any pip-audit
# vulnerability, and any grep-audit hit block. Accepted findings must be
# recorded in scripts/package/security_baseline.txt with an owner column and a
# dated comment explaining the acceptance. Low-severity bandit findings are
# reported for review but never block.
#
# All reports land under output/ for local inspection and CI artifacts.
#
# Baseline format (scripts/package/security_baseline.txt):
#   SEVERITY|OWNER|ID|LOCATION
#     SEVERITY: low | medium | high | vuln | pattern
#     OWNER:    responsible engineer or team (never omitted)
#     ID:       bandit test id (e.g. B108), PYSEC-/GHSA- id (pip-audit),
#               PATTERN:<name> for grep-audit hits, or `-` (any)
#     LOCATION: path substring, package name (pip-audit), or `-` (any)
# Lines starting with '#' are comments; every entry MUST carry a dated comment
# above it explaining why the finding is accepted.
set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

check_uv || exit 1

REPO_ROOT="$(get_repo_root)"
cd "$REPO_ROOT"

SRC_DIR="src/metainformant"
BASELINE_FILE="$SCRIPT_DIR/security_baseline.txt"
OUTPUT_DIR="output"
BANDIT_JSON="$OUTPUT_DIR/security_bandit.json"
PIP_AUDIT_JSON="$OUTPUT_DIR/security_pip_audit.json"
PATTERN_AUDIT_TXT="$OUTPUT_DIR/security_pattern_audit.txt"
PATTERN_TRAVERSAL_TXT="$OUTPUT_DIR/security_pattern_traversal.txt"
TMPDIR_FALLBACK="${TMPDIR:-/tmp}"
PIP_AUDIT_REQS="$(mktemp "$TMPDIR_FALLBACK/metainformant_requirements.XXXXXX.txt")"
trap 'rm -f "$PIP_AUDIT_REQS"' EXIT

mkdir -p "$OUTPUT_DIR"

run_bandit() {
    print_status "INFO" "Bandit: scanning $SRC_DIR (tests excluded)..."
    local rc=0
    uvx bandit -q -r "$SRC_DIR" -f json -o "$BANDIT_JSON" || rc=$?
    if [ "$rc" -eq 0 ]; then
        print_status "OK" "Bandit reported no findings"
    else
        print_status "INFO" "Bandit exited $rc; findings recorded in $BANDIT_JSON (verdict evaluated below)"
    fi
}

run_deps() {
    print_status "INFO" "pip-audit: exporting uv.lock-derived requirements (read-only)..."
    uv export --frozen --no-emit-project --format requirements-txt --no-hashes -o "$PIP_AUDIT_REQS"
    print_status "INFO" "pip-audit: scanning exported lockfile requirements..."
    local rc=0
    uvx pip-audit -r "$PIP_AUDIT_REQS" -f json -o "$PIP_AUDIT_JSON" || rc=$?
    if [ "$rc" -eq 0 ]; then
        print_status "OK" "pip-audit reported no known vulnerabilities"
    else
        print_status "INFO" "pip-audit exited $rc; vulnerabilities recorded in $PIP_AUDIT_JSON (verdict evaluated below)"
    fi
}

run_patterns() {
    print_status "INFO" "Pattern audit: shell-execution and path-traversal greps over $SRC_DIR..."
    local rc=0
    # Blocking shell-execution patterns: any hit means arguments or code reach
    # a shell/evaluator instead of an argv list. Pure comment lines cannot
    # execute and are filtered out.
    grep -rnE \
        --include='*.py' \
        --exclude-dir=__pycache__ \
        --exclude-dir=tests \
        -e 'shell[[:space:]]*=[[:space:]]*True' \
        -e 'os\.system\(' \
        -e 'os\.popen\(' \
        -e '(^|[^A-Za-z0-9_.])eval\(' \
        -e '(^|[^A-Za-z0-9_.])exec\(' \
        "$SRC_DIR" 2>/dev/null > "$PATTERN_AUDIT_TXT.raw" || rc=$?
    if [ "$rc" -ne 0 ] && [ "$rc" -ne 1 ]; then
        print_status "ERROR" "Pattern audit grep failed with exit $rc"
        exit "$rc"
    fi
    grep -v '^[^:]*:[0-9]*:[[:space:]]*#' "$PATTERN_AUDIT_TXT.raw" > "$PATTERN_AUDIT_TXT" || true
    rm -f "$PATTERN_AUDIT_TXT.raw"

    # Informational path-traversal candidates: quoted string literals that
    # begin with a parent-directory component. Benign occurrences include
    # sanitizer pattern lists and HTML templates with relative links; hits
    # are reported for review but never block.
    grep -rnE \
        --include='*.py' \
        --exclude-dir=__pycache__ \
        --exclude-dir=tests \
        -e "['\"]\.\./" \
        "$SRC_DIR" > "$PATTERN_TRAVERSAL_TXT" 2>/dev/null || true

    if [ -s "$PATTERN_AUDIT_TXT" ]; then
        print_status "WARN" "Pattern audit found $(wc -l < "$PATTERN_AUDIT_TXT" | tr -d ' ') shell-execution hit(s); verdict evaluated below"
    else
        print_status "OK" "Pattern audit found no shell-execution hits"
    fi
}

# Evaluate all produced reports against the baseline. Exit non-zero when any
# blocking finding is not covered by a baseline entry.
evaluate() {
    if [ ! -f "$BASELINE_FILE" ]; then
        print_status "ERROR" "Baseline file missing: $BASELINE_FILE"
        exit 1
    fi

    python3 - "$BANDIT_JSON" "$PIP_AUDIT_JSON" "$PATTERN_AUDIT_TXT" "$PATTERN_TRAVERSAL_TXT" "$BASELINE_FILE" <<'PY'
import json
import sys

bandit_path, pip_audit_path, pattern_path, traversal_path, baseline_path = sys.argv[1:6]

blocking = []
info_counts = {"bandit_low": 0}

baseline = []
with open(baseline_path, encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = [f.strip() for f in line.split("|")]
        if len(fields) != 4:
            sys.exit(f"Malformed baseline line (need SEVERITY|OWNER|ID|LOCATION): {line!r}")
        baseline.append(fields)


def baselined(severity, finding_id, location):
    for sev, _owner, bid, bloc in baseline:
        if sev not in ("*", severity):
            continue
        if bid not in ("-", "*", finding_id):
            continue
        if bloc not in ("-", "*") and bloc not in location:
            continue
        return True
    return False


try:
    with open(bandit_path, encoding="utf-8") as fh:
        bandit = json.load(fh)
    for result in bandit.get("results", []):
        severity = result["issue_severity"].lower()
        finding_id = result["test_id"]
        location = f"{result['filename']}:{result['line_number']}"
        if severity in ("medium", "high"):
            if not baselined(severity, finding_id, result["filename"]):
                blocking.append(
                    f"bandit {severity.upper()} {finding_id} "
                    f"{result['test_name']} at {location}"
                )
        else:
            info_counts["bandit_low"] += 1
except FileNotFoundError:
    print("INFO: bandit report not produced (skipped)")

try:
    with open(pip_audit_path, encoding="utf-8") as fh:
        pip_audit = json.load(fh)
    for dep in pip_audit.get("dependencies", []):
        for vuln in dep.get("vulns", []):
            vuln_id = vuln.get("id", "unknown")
            if not baselined("vuln", vuln_id, dep["name"]):
                fixes = ", ".join(vuln.get("fix_versions", []) or ["none"])
                blocking.append(f"pip-audit vuln {vuln_id} on {dep['name']}=={dep['version']} (fix: {fixes})")
except FileNotFoundError:
    print("INFO: pip-audit report not produced (skipped)")

try:
    with open(pattern_path, encoding="utf-8") as fh:
        for line in fh:
            hit = line.rstrip("\n")
            if not hit:
                continue
            path = hit.split(":", 1)[0]
            if not baselined("pattern", "PATTERN", path):
                blocking.append(f"pattern audit hit: {hit}")
except FileNotFoundError:
    print("INFO: pattern audit report not produced (skipped)")

try:
    with open(traversal_path, encoding="utf-8") as fh:
        candidates = [line.rstrip("\n") for line in fh if line.strip()]
    for hit in candidates:
        print(f"REVIEW (non-blocking): path-traversal candidate: {hit}")
except FileNotFoundError:
    candidates = []
    print("INFO: path-traversal report not produced (skipped)")

print(
    f"Baseline entries: {len(baseline)}; "
    f"bandit low-severity (non-blocking): {info_counts['bandit_low']}; "
    f"path-traversal review candidates (non-blocking): {len(candidates)}"
)
if blocking:
    print(f"BLOCKING: {len(blocking)} security finding(s) not covered by the baseline:")
    for finding in blocking:
        print(f"  - {finding}")
    sys.exit(1)
print("No blocking security findings.")
PY
}

case "${1:-all}" in
    "bandit")
        run_bandit
        evaluate
        ;;
    "deps")
        run_deps
        evaluate
        ;;
    "patterns")
        run_patterns
        evaluate
        ;;
    "all"|*)
        run_bandit
        run_deps
        run_patterns
        evaluate
        print_status "OK" "Security checks completed"
        ;;
esac
