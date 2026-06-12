#!/usr/bin/env bash
set -euo pipefail

echo "=== Spateo Environment Setup ==="
echo "Python: $(python3 --version 2>&1 || echo 'not found')"
echo "Platform: $(uname -s) $(uname -m)"

# Detect package manager
if command -v conda &>/dev/null; then
    MANAGER="conda"
    echo "Detected: conda"
elif command -v pip3 &>/dev/null || command -v pip &>/dev/null; then
    MANAGER="venv"
    echo "Detected: pip/venv (no conda)"
else
    echo "ERROR: No Python package manager found"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REQUIREMENTS="$SCRIPT_DIR/requirements.txt"

if [ "$MANAGER" = "conda" ]; then
    echo "Creating conda environment 'spateo'..."
    conda create -n spateo python=3.10 -y
    conda activate spateo
    echo "Installing requirements..."
    pip install --upgrade pip
    pip install -r "$REQUIREMENTS"
else
    VENV_DIR=".venv-spateo"
    echo "Creating venv at $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
    source "$VENV_DIR/bin/activate"
    echo "Installing requirements..."
    pip install --upgrade pip
    pip install -r "$REQUIREMENTS"
fi

echo ""
echo "=== Verification ==="
python -c "import spateo; print(f'spateo: {spateo.__version__}')" 2>&1 || echo "  spateo: FAILED"
python -c "import pyvista; print(f'pyvista: {pyvista.__version__}')" 2>&1 || echo "  pyvista: FAILED"
python -c "import scanpy; print(f'scanpy: {scanpy.__version__}')" 2>&1 || echo "  scanpy: FAILED"
python -c "import anndata; print(f'anndata: {anndata.__version__}')" 2>&1 || echo "  anndata: FAILED"

echo ""
echo "=== Setup Complete ==="
if [ "$MANAGER" = "conda" ]; then
    echo "Activate with: conda activate spateo"
else
    echo "Activate with: source $VENV_DIR/bin/activate"
fi
