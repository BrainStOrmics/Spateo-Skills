#!/usr/bin/env bash
# Sync Spateo-Skills to a Codex plugin directory.
# Usage: ./scripts/sync-to-codex-plugin.sh <target-dir>
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <target-plugin-dir>"
    exit 1
fi

TARGET="$1"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "=== Spateo-Skills → Codex Plugin Sync ==="
echo "Source: $REPO_ROOT"
echo "Target: $TARGET"

mkdir -p "$TARGET/skills" "$TARGET/docs"

rsync -av --exclude '__pycache__' --exclude '*.pyc' \
    "$REPO_ROOT/skills/" "$TARGET/skills/"
rsync -av "$REPO_ROOT/docs/" "$TARGET/docs/"

cp "$REPO_ROOT/README.md" "$TARGET/README.md"
cp "$REPO_ROOT/LICENSE" "$TARGET/LICENSE" 2>/dev/null || true
cp "$REPO_ROOT/pyproject.toml" "$TARGET/pyproject.toml" 2>/dev/null || true

cat > "$TARGET/plugin.json" << MANIFEST
{
  "name": "spateo-skills",
  "version": "1.0.0",
  "description": "Spateo spatial transcriptomics pipeline skills",
  "skills_directory": "skills",
  "minimum_cli_version": "2.0.0"
}
MANIFEST

echo "Sync complete. Files: $(find "$TARGET" -type f | wc -l)"
