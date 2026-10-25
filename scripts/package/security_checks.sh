#!/bin/bash
# Blocking security and dependency checks for the MetaInformAnt release surface.
#
# Mirrors the blocking security-checks job in .github/workflows/security.yml so
# developers can reproduce the CI gate locally before pushing. Three gates:
#
#   1. bandit static analysis over src/metainformant (tests excluded)
#   2. pip-audit over every package pinned in uv.lock: all extras, dependency
#      groups and environment markers, exported as a PEP 751 pylock.toml and
#      audited in place against OSV, which carries both GitHub (GHSA) and PyPA
#      (PYSEC) advisories (no pip, no virtual environment, nothing installed)
#   3. grep audit for shell-execution and path-traversal patterns in
#      maintained src/ paths
#
# Gate policy: medium-and-higher severity bandit findings, any pip-audit
# vulnerability, any locked dependency pip-audit could not audit, and any
# grep-audit hit block. A scanner that runs but leaves no complete report also
# blocks: a failed scan is missing evidence, never a clean result. That covers
# bandit's own read/parse errors and a zero-line scan, and a pip-audit report
# that does not account for every package in the exported lock. Accepted findings
# must be recorded in scripts/package/security_baseline.txt with an owner
# column and a dated comment explaining the acceptance. Low-severity bandit
# findings are reported for review but never block.
#
# All reports land under output/ for local inspection and CI artifacts. Each
# scanner removes its previous report first, and only the scanners run by this
# invocation are evaluated, so a stale report can never stand in for a scan.
#
# Baseline format (scripts/package/security_baseline.txt):
#   SEVERITY|OWNER|ID|LOCATION
#     SEVERITY: low | medium | high | vuln | pattern
#     OWNER:    responsible engineer or team (never omitted)
#     ID:       bandit test id (e.g. B108), PYSEC-/GHSA- id (pip-audit),
#               `unaudited` for a locked dependency pip-audit skipped,
#               PATTERN for grep-audit hits, or `-` (any)
#     LOCATION: path substring (bandit, pattern), exact package name
#               (pip-audit), or `-` (any)
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
PIP_AUDIT_LOCK_DIR="$(mktemp -d "$TMPDIR_FALLBACK/metainformant_pylock.XXXXXX")"
trap 'rm -rf "$PIP_AUDIT_LOCK_DIR"' EXIT

mkdir -p "$OUTPUT_DIR"

run_bandit() {
    print_status "INFO" "Bandit: scanning $SRC_DIR (tests excluded)..."
    rm -f "$BANDIT_JSON"
    local rc=0
    uvx bandit -q -r "$SRC_DIR" -f json -o "$BANDIT_JSON" || rc=$?
    if [ "$rc" -eq 0 ]; then
        print_status "OK" "Bandit reported no findings"
    else
        print_status "INFO" "Bandit exited $rc; verdict evaluated below from $BANDIT_JSON"
    fi
}

run_deps() {
    print_status "INFO" "pip-audit: exporting every uv.lock package to pylock.toml (read-only)..."
    rm -f "$PIP_AUDIT_JSON"
    # uv also echoes the exported lock on stdout; the written file is the artifact.
    uv export --frozen --no-emit-project --all-extras --all-groups --format pylock.toml \
        -o "$PIP_AUDIT_LOCK_DIR/pylock.toml" > /dev/null
    print_status "INFO" "pip-audit: auditing the exported lock in place..."
    local rc=0
    uvx --from 'pip-audit>=2.10' pip-audit --locked "$PIP_AUDIT_LOCK_DIR" -s osv -f json -o "$PIP_AUDIT_JSON" || rc=$?
    if [ "$rc" -eq 0 ]; then
        print_status "OK" "pip-audit reported no known vulnerabilities"
    else
        print_status "INFO" "pip-audit exited $rc; verdict evaluated below from $PIP_AUDIT_JSON"
    fi
}

run_patterns() {
    print_status "INFO" "Pattern audit: shell-execution and path-traversal greps over $SRC_DIR..."
    rm -f "$PATTERN_AUDIT_TXT" "$PATTERN_TRAVERSAL_TXT"
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
    local filter_rc=0
    grep -v '^[^:]*:[0-9]*:[[:space:]]*#' "$PATTERN_AUDIT_TXT.raw" > "$PATTERN_AUDIT_TXT" || filter_rc=$?
    if [ "$filter_rc" -gt 1 ]; then
        print_status "ERROR" "Pattern audit comment filter failed with exit $filter_rc"
        exit "$filter_rc"
    fi
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

# Evaluate the reports of the named scanners (those run by this invocation)
# against the baseline with security_evaluate.py. Exit non-zero when a named
# scanner left no complete report or a blocking finding is not baselined.
evaluate() {
    if [ ! -f "$BASELINE_FILE" ]; then
        print_status "ERROR" "Baseline file missing: $BASELINE_FILE"
        exit 1
    fi

    local scanners
    scanners=$(IFS=,; echo "$*")
    python3 "$SCRIPT_DIR/security_evaluate.py" \
        --ran "$scanners" \
        --baseline "$BASELINE_FILE" \
        --bandit "$BANDIT_JSON" \
        --pip-audit "$PIP_AUDIT_JSON" \
        --pylock "$PIP_AUDIT_LOCK_DIR/pylock.toml" \
        --patterns "$PATTERN_AUDIT_TXT" \
        --traversal "$PATTERN_TRAVERSAL_TXT"
}

case "${1:-all}" in
    "bandit")
        run_bandit
        evaluate bandit
        ;;
    "deps")
        run_deps
        evaluate deps
        ;;
    "patterns")
        run_patterns
        evaluate patterns
        ;;
    "all"|*)
        run_bandit
        run_deps
        run_patterns
        evaluate bandit deps patterns
        print_status "OK" "Security checks completed"
        ;;
esac
