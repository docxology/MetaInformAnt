#!/bin/bash
# Generalized script to set up temp directories for any repo location
# Works on external drives, home directories, or any filesystem
# Location: scripts/core/setup_temp_dirs.sh

set -e  # Exit on error

# Determine repository root (works from any subdirectory)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "🔧 Setting up repository-local temp directories..."
echo "   Repository: $REPO_ROOT"

# Detect filesystem type
FS_TYPE=$(stat -f -c %T "$REPO_ROOT" 2>/dev/null || echo "unknown")
echo "   Filesystem: $FS_TYPE"

# Create temp directory structure
TEMP_ROOT="$REPO_ROOT/.tmp"
mkdir -p "$TEMP_ROOT/bash"
mkdir -p "$TEMP_ROOT/python"
mkdir -p "$TEMP_ROOT/git"
mkdir -p "$TEMP_ROOT/downloads"

# Set permissions (like /tmp) - works on all filesystems
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
echo "ℹ️  System /tmp status:"
df -h /tmp 2>/dev/null | tail -1 || echo "   (Unable to check /tmp)"

# Test bash heredoc (commonly fails when /tmp is full)
echo ""
echo "🧪 Testing bash heredoc..."
if cat <<EOF > /dev/null 2>&1
Test heredoc
EOF
then
    echo "✅ Bash heredoc works!"
else
    echo "❌ Bash heredoc still fails"
    exit 1
fi

# Generate .bashrc snippet (generalized)
echo ""
echo "📝 To make this permanent, add to ~/.bashrc:"
echo ""
cat << 'BASHRC_SNIPPET'
# METAINFORMANT repository-local temp directories
# Adjust REPO_PATH to match your repository location
REPO_PATH="$HOME/MetaInformAnt"  # Change this to your repo path

if [ -d "$REPO_PATH" ]; then
    export TMPDIR="$REPO_PATH/.tmp/bash"
    export TEMP="$TMPDIR"
    export TMP="$TMPDIR"
    export UV_CACHE_DIR="$REPO_PATH/.cache/uv"
    export PIP_CACHE_DIR="$REPO_PATH/.cache/pip"
fi
BASHRC_SNIPPET

echo ""
echo "Or for this specific location:"
echo ""
echo "# METAINFORMANT temp directories"
echo "if [ -d \"$REPO_ROOT\" ]; then"
echo "    export TMPDIR=\"$TEMP_ROOT/bash\""
echo "    export TEMP=\"\$TMPDIR\""
echo "    export TMP=\"\$TMPDIR\""
echo "    export UV_CACHE_DIR=\"$CACHE_ROOT/uv\""
echo "    export PIP_CACHE_DIR=\"$CACHE_ROOT/pip\""
echo "fi"
echo ""
echo "✅ Setup complete! You can now run commands without temp directory issues."
echo ""
echo "💡 Tip: Run 'source ~/.bashrc' after adding to .bashrc to apply immediately."

