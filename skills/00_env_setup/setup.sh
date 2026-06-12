#!/usr/bin/env bash
set -euo pipefail

# ── Spateo Environment Setup ──────────────────────────────────────
# Detects conda and uv, prompts user to choose, then installs.
# ───────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV="environment.yml"
UV_REQUIREMENTS="uv_requirements.txt"
UV_PYPROJECT="uv_pyproject.toml"
ENV_NAME="spateo"

echo "========================================"
echo "  Spateo Environment Setup"
echo "========================================"
echo ""

# Detect available tools
HAS_CONDA=false
HAS_UV=false
PYTHON_VER=""

if command -v conda &>/dev/null; then
    HAS_CONDA=true
    echo "[OK] conda: $(conda --version 2>&1 | head -1)"
fi

if command -v uv &>/dev/null; then
    HAS_UV=true
    echo "[OK] uv:     $(uv --version 2>&1)"
fi

if python3 --version &>/dev/null; then
    PYTHON_VER="$(python3 --version 2>&1)"
    echo "[OK] python: $PYTHON_VER"
fi
echo ""

# If neither, abort
if [ "$HAS_CONDA" = false ] && [ "$HAS_UV" = false ]; then
    echo "ERROR: Neither conda nor uv found."
    echo "Install one of:"
    echo "  conda: https://docs.conda.io/en/latest/miniconda.html"
    echo "  uv:    https://docs.astral.sh/uv/"
    exit 1
fi

# Auto-select if only one available
if [ "$HAS_CONDA" = true ] && [ "$HAS_UV" = false ]; then
    CHOICE="conda"
elif [ "$HAS_CONDA" = false ] && [ "$HAS_UV" = true ]; then
    CHOICE="uv"
else
    # Both available — ask user
    echo "Both conda and uv detected. Choose environment manager:"
    echo "  1) conda  (full dependency resolution, larger env)"
    echo "  2) uv     (fast, lightweight, pip-compatible)"
    echo ""
    read -rp "Select [1/2] (default: 2): " ANSWER
    case "$ANSWER" in
        1) CHOICE="conda" ;;
        *) CHOICE="uv" ;;
    esac
fi

echo ""
echo "Using: $CHOICE"
echo "========================================"
echo ""

if [ "$CHOICE" = "conda" ]; then
    echo "Creating conda environment '$ENV_NAME'..."
    conda env create -f "$SCRIPT_DIR/$CONDA_ENV" -y
    echo ""
    echo "========================================"
    echo "  Setup Complete"
    echo "========================================"
    echo "Activate with: conda activate $ENV_NAME"
    echo ""
    echo "Verify:"
    echo "  conda activate $ENV_NAME"
    echo "  python -c 'import spateo; print(spateo.__version__)'"

elif [ "$CHOICE" = "uv" ]; then
    echo "Creating uv-managed venv..."

    # Use uv venv to create environment
    VENV_DIR=".venv-spateo"
    if command -v python3.10 &>/dev/null; then
        uv venv "$VENV_DIR" --python 3.10
    else
        uv venv "$VENV_DIR" --python 3.10
    fi

    # Activate venv
    source "$VENV_DIR/bin/activate"

    # Install via uv pip (fast, uses requirements.txt)
    echo "Installing dependencies (uv pip)..."
    uv pip install -r "$SCRIPT_DIR/$UV_REQUIREMENTS"

    echo ""
    echo "========================================"
    echo "  Setup Complete"
    echo "========================================"
    echo "Activate with: source $VENV_DIR/bin/activate"
    echo ""
    echo "Verify:"
    echo "  source $VENV_DIR/bin/activate"
    echo "  python -c 'import spateo; print(spateo.__version__)'"
    echo ""
    echo "Alternative: use pyproject.toml"
    echo "  uv sync --project $SCRIPT_DIR/$UV_PYPROJECT"
fi
