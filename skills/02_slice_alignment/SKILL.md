---
name: spateo-slice-alignment
description: Spatial alignment of tissue slices — 2D→3D registration
---

# Stage 1: Slice Alignment

Register multiple tissue slices into a common coordinate space using `morpho_align` and `morpho_align_ref`.

## Quick Start

```bash
# Align two slices
python -m skills.02_slice_alignment.two_slice \
  --slice1 ./data/slice_a.h5ad --slice2 ./data/slice_b.h5ad

# Align multiple slices (z-splitting + 3D concat)
python -m skills.02_slice_alignment.multi_slice \
  --input ./data/slices_dir/ --outdir ./output/

# Two-stage 3D alignment (point cloud + reference-based)
python -m skills.02_slice_alignment.two_stage_3d \
  --stage1 ./data/stage1.h5ad --stage2 ./data/stage2.h5ad
```

## Modules

### `two_slice.py` — Two-slice alignment

```python
from skills.02_slice_alignment.two_slice import (
    TwoSliceConfig, TwoSliceResult, run_two_slice_pipeline,
)

config = TwoSliceConfig(
    slice1_path="./data/slice_a.h5ad",
    slice2_path="./data/slice_b.h5ad",
    preprocess=True,
    use_downsampling=False,
    max_iter=300,
)
result = run_two_slice_pipeline(config)
```

**Key APIs:** `st.align.morpho_align`, `st.align.morpho_align_ref`

### `multi_slice.py` — Multi-slice continuous alignment

Splits data by z-coordinate, aligns consecutive pairs, concatenates in 3D.

```python
from skills.02_slice_alignment.multi_slice import (
    MultiSliceConfig, run_multi_slice_pipeline,
)

config = MultiSliceConfig(
    input_path="./data/slices/",
    z_key="z_coord",
    groupby="tissue_type",
    max_iter=300,
)
result = run_multi_slice_pipeline(config)
```

### `two_stage_3d.py` — Two-stage 3D reconstruction + alignment

Stage 1: Build point cloud/mesh from individual slices.
Stage 2: Align using reference-based `morpho_align_ref`.

```python
from skills.02_slice_alignment.two_stage_3d import (
    TwoStage3DConfig, run_two_stage_3d_pipeline,
)
```

## Input/Output

**Input:** AnnData files with `.obsm["spatial"]` (2D coordinates per slice).
**Output:** AnnData with `.obsm["align_spatial"]` (registered coordinates) or concatenated 3D AnnData.

## Key Parameters

| Parameter | Purpose | Default |
|-----------|---------|---------|
| `slice1_path`, `slice2_path` | Input slice paths | — |
| `preprocess` | Run spateo preprocessing pipeline | `True` |
| `spatial_key` | Key for input coordinates | `"spatial"` |
| `key_added` | Key for output coordinates | `"align_spatial"` |
| `max_iter` | Alignment iterations | 300 |
| `use_downsampling` | Use reference-based alignment | `False` |
| `beta` | Spatial smoothness weight | 1.0 |
| `lambdaVF` | Vector field regularization | 1.0 |
| `device` | CPU/GPU device | Auto-detected |

## Troubleshooting

| Problem | Fix |
|---------|-----|
| OOM on GPU | Set `use_downsampling=True` or `device="cpu"` |
| Missing spatial key | Verify `.obsm["spatial"]` exists in input |
| NaN in aligned coords | Check for NaN/Inf in input spatial coords |
