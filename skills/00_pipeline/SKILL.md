---
name: spateo-full-pipeline
description: Run the complete Spateo pipeline — IO through Morphogenesis
---

# Full Pipeline Orchestrator

Run all Spateo stages sequentially: `IO → Alignment → 3D Reconstruction → Morphogenesis`.

## Quick Start

```bash
# Run full pipeline with default settings
python -m skills.00_pipeline.run_all \
  --data-path ./data/xenium_outs/ \
  --platform xenium \
  --stage2-path ./data/stage2.h5ad \
  --out-dir ./output/

# Run with custom alignment
python -m skills.00_pipeline.run_all \
  --data-path ./data/ \
  --auto-detect \
  --slices-dir ./data/slices/ \
  --out-dir ./output/ \
  --skip-morphogenesis
```

## Pipeline Stages

```
Stage 0 (IO) ──→ Stage 1 (Alignment) ─→ Stage 2 (Reconstruction)
                                                    ↓
                                     interpolation    morphology
                                                    ↓
                                        Stage 3 (Morphogenesis)
                                        vectorfield ──→ feature
```

| Stage | What it does | Output |
|-------|-------------|--------|
| **0: IO** | Read raw spatial data from platform | `adata.h5ad` with `.obsm["spatial"]` |
| **1: Alignment** | Register slices into common coordinate space | `aligned.h5ad` with `.obsm["align_spatial"]` |
| **2: Reconstruction** | Build 3D models + interpolate expression | VTK models, interpolated data |
| **3: Morphogenesis** | Infer vector fields + compute features | Vector field layers, DEG tables |

## CLI Parameters

| Parameter | Purpose | Default |
|-----------|---------|---------|
| `--data-path` | Raw data directory (Stage 0) | — |
| `--platform` | Platform name (e.g. `xenium`, `merfish`) | `--auto-detect` to guess |
| `--slices-dir` | Directory with slice `.h5ad` files (Stage 1) | Uses IO output |
| `--stage2-path` | Second stage data for morphogenesis (Stage 3) | — |
| `--out-dir` | Output directory for all stages | `./output/` |
| `--skip-alignment` | Skip Stage 1 | — |
| `--skip-reconstruction` | Skip Stage 2 | — |
| `--skip-morphogenesis` | Skip Stage 3 | — |
| `--model-type` | 3D reconstruction model type | `surface` |
| `--interpolation-method` | Gene interpolation method | `vtk` |
| `--device` | Compute device | Auto-detect |
| `--verbose` | Enable debug logging | — |

## Programmatic Usage

```python
from skills.00_pipeline.run_all import PipelineConfig, run_full_pipeline

config = PipelineConfig(
    data_path="./data/xenium_outs/",
    platform="xenium",
    stage2_path="./data/stage2.h5ad",
    out_dir="./output/",
)
result = run_full_pipeline(config)
```

## Troubleshooting

| Problem | Fix |
|---------|-----|
| Stage 0 fails | Check data directory contains platform-specific files |
| Stage 1 fails | Verify `.obsm["spatial"]` exists in input AnnData |
| Stage 2 OOM | Reduce `voxel_size` or use `--model-type pointcloud` |
| Stage 3 fails | Ensure stage1 and stage2 have matching gene sets |
| Any stage hangs | Set `--device cpu` to avoid GPU issues |
