---
name: spateo-env-setup
description: Set up Spateo runtime environment — choose conda or uv
---

# Spateo Environment Setup

Configure a reproducible Python environment for Spateo spatial transcriptomics skills.
Supports both **conda** (full dependency resolution) and **uv** (fast, lightweight).

## Quick Setup

### Automated (interactive choice)

```bash
bash skills/00_env_setup/setup.sh
```

Script detects available tools and prompts you to choose between conda and uv.
If only one is installed, it auto-selects.

### Manual — conda

```bash
conda env create -f skills/00_env_setup/environment.yml
conda activate spateo
```

Best for:
- Complex native dependencies (geospatial, VTK, Open3D)
- GPU environments (pytorch-cuda included)
- Platforms where pip wheels fail (Linux ARM, macOS arm64)

### Manual — uv (recommended for speed)

```bash
# Option A: from requirements.txt (fastest)
uv venv .venv-spateo --python 3.10
source .venv-spateo/bin/activate
uv pip install -r skills/00_env_setup/uv_requirements.txt

# Option B: from pyproject.toml (locked resolution)
uv sync --project skills/00_env_setup/uv_pyproject.toml
```

Best for:
- Fast installation (10x faster than pip)
- CI/CD pipelines
- Developers who want a clean, reproducible venv

## Environment Files

| File | Manager | Description |
|------|---------|-------------|
| `environment.yml` | conda | Full conda-forge spec with pip fallback for spateo-release |
| `uv_requirements.txt` | uv/pip | Flat requirements, version ranges for Python 3.10 |
| `uv_pyproject.toml` | uv | pyproject.toml with `[tool.uv]` config and git sources |

## Key Dependencies

- **spateo-release** — core library (git install, not on PyPI)
- **torch** — deep learning backend (>=2.0,<2.5 for Python 3.10)
- **pyvista** — 3D visualization (requires VTK>=9.2)
- **scanpy** — single-cell analysis
- **anndata** — AnnData format
- **numpy** — pinned >=1.21,<1.24 for Spateo compatibility

## Verification

```bash
python -c "import spateo; print(f'spateo: {spateo.__version__}')"
python -c "import pyvista; print(f'pyvista: {pyvista.__version__}')"
python -c "import scanpy; print(f'scanpy: {scanpy.__version__}')"
```

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `numpy` version conflict | Install `numpy>=1.21,<1.24` first |
| `vtk` build failure | Use conda: `conda install -c conda-forge vtk` |
| `pyvista` headless | `export PYVISTA_OFF_SCREEN=true` |
| POT library missing | conda: `conda install -c conda-forge pot` |
| `dynamo-release` not on PyPI | conda: `conda install -c conda-forge dynamo-release` |
| uv too old | `uv self update` |
