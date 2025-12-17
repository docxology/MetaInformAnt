#!/bin/bash
# DEPRECATED: Use scripts/package/setup.sh instead
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
show_deprecation_warning "setup_uv.sh" "scripts/package/setup.sh"
exec "$SCRIPT_DIR/setup.sh" "$@"