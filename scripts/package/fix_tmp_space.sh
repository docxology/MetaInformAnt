#!/bin/bash
# Quick fix for /tmp space issues by cleaning pytest temp files and uv cache

echo "=== Checking /tmp usage ==="
df -h /tmp

echo ""
echo "=== Cleaning pytest temp directories ==="
if [ -d "/tmp/pytest-of-q" ]; then
    du -sh /tmp/pytest-of-q
    echo "Removing pytest temp files..."
    rm -rf /tmp/pytest-of-q
    echo "✓ Cleaned pytest temp files"
else
    echo "No pytest temp files found"
fi

echo ""
echo "=== Cleaning uv cache ==="
if [ -d "/tmp/uv-cache" ]; then
    du -sh /tmp/uv-cache
    echo "Removing uv cache..."
    rm -rf /tmp/uv-cache
    echo "✓ Cleaned uv cache"
else
    echo "No uv cache found"
fi

echo ""
echo "=== Finding other large temp files ==="
find /tmp -maxdepth 1 -type f -size +100M -exec ls -lh {} \; 2>/dev/null | head -10

echo ""
echo "=== After cleanup ==="
df -h /tmp

