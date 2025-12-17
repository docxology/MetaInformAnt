#!/bin/bash
# DEPRECATED: Use scripts/package/verify.sh --mode setup instead
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"
show_deprecation_warning "verify_uv_setup.sh" "scripts/package/verify.sh --mode setup"
exec "$SCRIPT_DIR/verify.sh" --mode setup "$@"