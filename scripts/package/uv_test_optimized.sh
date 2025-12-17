#!/bin/bash
# DEPRECATED: Use scripts/package/test.sh --mode fast instead
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
show_deprecation_warning "uv_test_optimized.sh" "scripts/package/test.sh --mode fast"
exec "$SCRIPT_DIR/test.sh" --mode fast "$@"