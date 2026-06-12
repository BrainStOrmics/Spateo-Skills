---
name: spateo-data-io
description: Read spatial transcriptomics data — 10 platforms + auto-detect
---

# Stage 0: Data IO

Load spatial transcriptomics data from 10+ platforms into standardized `AnnData` format with `.obsm["spatial"]`.

## Quick Start

```bash
# Manual read — specify platform
python -m skills.01_data_io.manual_read --platform xenium --data-path ./data/xenium_outs/

# Auto-detect — let the tool figure out the platform
python -m skills.01_data_io.auto_read --data-path ./data/unknown_dataset/
```

## Modules

### `manual_read.py` — Platform-specific readers

Dispatch by platform name. Each reader wraps `spateo.io.protocol_io.spatial.*`.

```python
from skills.01_data_io.manual_read import read_by_platform

adata = read_by_platform("xenium", "./data/xenium_outs/")
adata = read_by_platform("merfish", "./data/merfish/", counts_file="cell_by_gene.csv")
```

**Supported platforms:**

| Key | Platform | Key File(s) |
|-----|----------|-------------|
| `merfish` | MERFISH | `cell_by_gene.csv`, `cell_metadata.csv` |
| `seqfish` | seqFISH | CxG CSV, coordinates CSV, DAPI TIFF |
| `slideseq` | Slide-seq | Standard output directory |
| `starmap_plus` | STARmap+ | Expression/spatial CSVs |
| `stereoseq` | Stereo-seq | GEM file (binsize default=50) |
| `stereoseq_agg` | Stereo-seq agg | GEM file + image file |
| `xenium` | Xenium | 10x outs directory |
| `visium` | Visium | 10x Visium output |
| `visium_hd` | Visium HD | 10x Visium HD binned output |
| `openst` | Open-ST | Open-ST output |

### `auto_read.py` — Auto-detect + batch read

```python
from skills.01_data_io.auto_read import (
    detect_spatial_platform,
    read_auto_spatial_data,
    read_many_auto,
)

# Detect platform type
platform = detect_spatial_platform("./data/unknown/")

# Auto-detect + read
adata = read_auto_spatial_data("./data/unknown/")

# Batch read multiple directories
results = read_many_auto(["./data/dataset1/", "./data/dataset2/"])
```

## Input/Output

**Input:** Platform-specific raw data directories/files.
**Output:** `AnnData` object with `.obsm["spatial"]` containing 2D/3D coordinates.

## Key Parameters

| Parameter | Purpose | Default |
|-----------|---------|---------|
| `data_path` | Path to raw data directory | — |
| `platform` | Platform name (manual_read) | — |
| `counts_file` | Gene expression matrix filename | `cell_by_gene.csv` |
| `binsize` | Stereo-seq binning size | 50 |
| `load_images` | Load DAPI/tissue images | `True` |

## Troubleshooting

| Problem | Fix |
|---------|-----|
| ImportError: spateo | Run `bash skills/00_env_setup/setup.sh` |
| File not found | Check `data_path` is absolute or relative to working dir |
| Platform not detected | Verify input directory contains expected files for the platform |
