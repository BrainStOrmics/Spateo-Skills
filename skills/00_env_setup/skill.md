---
name: spateo-env-setup
description: Set up Spateo runtime environment with conda or venv
---

# Spateo Environment Setup

## Purpose
Configure a reproducible Python environment for running Spateo spatial transcriptomics skills.

## Prerequisites
- Python 3.10 (conda or system)
- Linux (x86_64) or macOS (arm64/x86_64)

## Method 1: Conda (recommended for complex deps)

```bash
conda create -n spateo python=3.10 -y
conda activate spateo
pip install -r Spateo-Skills/00_env_setup/requirements.txt
```

## Method 2: venv + pip

```bash
python3.10 -m venv .venv-spateo
source .venv-spateo/bin/activate
pip install --upgrade pip
pip install -r Spateo-Skills/00_env_setup/requirements.txt
```

## Key Dependencies
- **spateo-release** — core library (git install)
- **torch** — deep learning backend
- **pyvista** — 3D visualization (requires VTK)
- **scanpy** — single-cell analysis
- **anndata** — AnnData format
- **numpy==1.23.5** — pinned for spateo compatibility
- **scipy==1.10.1** — scientific computing

## Platform Notes
- **Linux**: VTK/pyvista install cleanly via pip
- **macOS arm64**: Some packages (open3d, vtk) may need conda-forge
- The original `envs/environment.yml` was macOS-specific with hardcoded prefix
- `requirements.txt` in this directory is cleaned for cross-platform use

## Verification

```bash
python -c "import spateo; print(spateo.__version__)"
python -c "import pyvista; print(pyvista.__version__)"
python -c "import scanpy; print(scanpy.__version__)"
```

## Troubleshooting
- `numpy` version conflict: install `numpy==1.23.5` first, then other deps
- `vtk` build failure: install `vtk` via conda (`conda install -c conda-forge vtk`)
- `pyvista` headless: `export PYVISTA_OFF_SCREEN=true`
- POT library: `conda install -c conda-forge pot` if pip fails
