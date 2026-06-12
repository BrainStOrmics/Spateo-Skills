---
name: spateo-3d-reconstruction
description: 3D model generation, gene expression interpolation, morphology analysis
---

# Stage 2: 3D Reconstruction

Generate 3D models (point cloud, surface, mesh, voxel) from aligned slices, interpolate gene expression onto 3D grids, and compute morphology metrics.

## Quick Start

```bash
# 3D model reconstruction
python -m skills.03_3d_reconstruction.reconstruction \
  --input ./data/aligned.h5ad --out-dir ./output --groupby tissue_type

# Gene expression interpolation onto 3D voxel grid
python -m skills.03_3d_reconstruction.interpolation \
  --input ./data/aligned.h5ad --method vtk --out-dir ./output

# Morphology analysis (KDE density + metrics)
python -m skills.03_3d_reconstruction.morphology \
  --input ./data/aligned.h5ad --out-dir ./output
```

## Modules

### `reconstruction.py` — 3D model generation

Creates VTK/PyVista models: point cloud, surface mesh, cell mesh, voxel grid, subtype segmentation.

```python
from skills.03_3d_reconstruction.reconstruction import (
    ReconstructionConfig, run_reconstruction_pipeline,
)

config = ReconstructionConfig(
    input_path="./data/aligned.h5ad",
    model_type="surface",
    groupby="tissue_type",
    out_dir="./output/",
)
result = run_reconstruction_pipeline(config)
```

**Model types:** `pointcloud`, `surface`, `cell_mesh`, `voxel`, `subtype`

### `interpolation.py` — Gene expression interpolation

Interpolates sparse gene expression onto dense 3D voxel grids.

```python
from skills.03_3d_reconstruction.interpolation import (
    InterpolationConfig, run_interpolation_pipeline,
)

config = InterpolationConfig(
    input_path="./data/aligned.h5ad",
    method="vtk",  # vtk | gp | kernel | deep
    voxel_size=10.0,
    gene_list=["Sox2", "Pax6"],
)
```

**Methods:**

| Method | Description | Speed | Accuracy |
|--------|-------------|-------|----------|
| `vtk` | VTK-based interpolation | Fast | Good |
| `gp` | Gaussian Process | Slow | Best |
| `kernel` | Kernel density | Medium | Good |
| `deep` | Deep learning | Slow (GPU) | Best |

### `morphology.py` — Morphology metrics

Computes KDE-based spatial density and morphology metrics.

```python
from skills.03_3d_reconstruction.morphology import (
    MorphologyConfig, run_morphology_pipeline,
)
```

## Input/Output

**Input:** Aligned AnnData with `.obsm["align_spatial"]` or 3D spatial coords.
**Output:** VTK model files (`.vtk`, `.ply`, `.stl`), interpolated expression arrays, morphology metric tables.

## Key Parameters

| Parameter | Purpose | Default |
|-----------|---------|---------|
| `input_path` | Aligned AnnData path | — |
| `model_type` | 3D model type | `"surface"` |
| `method` | Interpolation method | `"vtk"` |
| `groupby` | Column for grouping cells | `"tissue_type"` |
| `voxel_size` | Voxel grid resolution | 10.0 |
| `gene_list` | Genes to interpolate | All genes |

## Troubleshooting

| Problem | Fix |
|---------|-----|
| VTK import error | Check pyvista+vtk installed |
| Headless display | `export PYVISTA_OFF_SCREEN=true` |
| OOM on interpolation | Reduce `voxel_size` or use `method="vtk"` |
| No spatial key found | Run alignment first (Stage 1) |
