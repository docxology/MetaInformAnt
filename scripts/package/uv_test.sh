#!/bin/bash
# DEPRECATED: Use scripts/package/test.sh instead
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
show_deprecation_warning "uv_test.sh" "scripts/package/test.sh"
exec "$SCRIPT_DIR/test.sh" "$@"