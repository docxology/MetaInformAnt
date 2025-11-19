#!/bin/bash
# Fix "out of space" errors by setting up external drive temp directories
# Location: scripts/core/fix_disk_space.sh

set -e  # Exit on error

# Determine repository root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "🔧 Setting up external drive temp directories..."
echo "   Repository: $REPO_ROOT"

# Create temp directory structure
TEMP_ROOT="$REPO_ROOT/.tmp"
mkdir -p "$TEMP_ROOT/bash"
mkdir -p "$TEMP_ROOT/python"
mkdir -p "$TEMP_ROOT/git"
mkdir -p "$TEMP_ROOT/downloads"

# Set permissions (like /tmp)
chmod 1777 "$TEMP_ROOT/bash" 2>/dev/null || chmod 777 "$TEMP_ROOT/bash"

echo "✅ Created temp directories:"
echo "   $TEMP_ROOT/bash"
echo "   $TEMP_ROOT/python"
echo "   $TEMP_ROOT/git"
echo "   $TEMP_ROOT/downloads"

# Create cache directory structure
CACHE_ROOT="$REPO_ROOT/.cache"
mkdir -p "$CACHE_ROOT/uv"
mkdir -p "$CACHE_ROOT/pip"
mkdir -p "$CACHE_ROOT/npm"
mkdir -p "$CACHE_ROOT/cargo"

echo "✅ Created cache directories:"
echo "   $CACHE_ROOT/uv"
echo "   $CACHE_ROOT/pip"
echo "   $CACHE_ROOT/npm"
echo "   $CACHE_ROOT/cargo"

# Export environment variables
export TMPDIR="$TEMP_ROOT/bash"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

echo ""
echo "✅ Environment variables set for current session:"
echo "   TMPDIR=$TMPDIR"
echo "   TEMP=$TEMP"
echo "   TMP=$TMP"

# Check current disk usage
echo ""
echo "📊 Disk usage:"
df -h "$REPO_ROOT" | head -2

# Check system /tmp status
echo ""
echo "⚠️  System /tmp status:"
df -h /tmp | tail -1

# Test bash heredoc (previously failing)
echo ""
echo "🧪 Testing bash heredoc (was failing before)..."
if cat <<EOF > /dev/null 2>&1
Test heredoc
EOF
then
    echo "✅ Bash heredoc works!"
else
    echo "❌ Bash heredoc still fails"
    exit 1
fi

# Instructions for persistence
echo ""
echo "📝 To make this permanent, add to ~/.bashrc:"
echo ""
echo "# METAINFORMANT external drive temp directories"
echo "if [ -d \"$REPO_ROOT\" ]; then"
echo "    export TMPDIR=\"$TEMP_ROOT/bash\""
echo "    export TEMP=\"\$TMPDIR\""
echo "    export TMP=\"\$TMPDIR\""
echo "    export UV_CACHE_DIR=\"$CACHE_ROOT/uv\""
echo "    export PIP_CACHE_DIR=\"$CACHE_ROOT/pip\""
echo "fi"
echo ""
echo "✅ Setup complete! You can now run commands without 'out of space' errors."

