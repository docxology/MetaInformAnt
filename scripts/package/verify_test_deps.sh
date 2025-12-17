#!/bin/bash
# DEPRECATED: Use scripts/package/verify.sh --mode deps instead
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
show_deprecation_warning "verify_test_deps.sh" "scripts/package/verify.sh --mode deps"
exec "$SCRIPT_DIR/verify.sh" --mode deps "$@"