# Spateo-Skills

Standardized Python skill modules for the **Spateo** spatial transcriptomics
pipeline — from raw data ingestion through 3D reconstruction to morphogenesis
analysis.

## Quickstart

### 1. Install Environment

```bash
# Option A: automated setup (prompts conda/uv choice if both available)
bash skills/00_env_setup/setup.sh

# Option B: conda
conda env create -f skills/00_env_setup/environment.yml
conda activate spateo

# Option C: uv (faster)
uv venv .venv-spateo --python 3.10 && source .venv-spateo/bin/activate
uv pip install -r skills/00_env_setup/uv_requirements.txt
```

### 2. Run the Full Pipeline

```bash
python -m skills.00_pipeline.run_all \
  --data-path ./data/xenium_outs/ \
  --platform xenium \
  --stage2-path ./data/stage2.h5ad \
  --out-dir ./output/
```

Or run stages individually:

### 2b. Run Individual Skills

# Stage 1 — Align two tissue slices
python -m skills.02_slice_alignment.two_slice \
  --slice1 ./data/slice_a.h5ad --slice2 ./data/slice_b.h5ad

# Stage 2 — 3D reconstruction
python -m skills.03_3d_reconstruction.reconstruction \
  --input ./data/aligned.h5ad --out-dir ./output --groupby tissue_type

# Stage 3 — Morphogenesis vector field
python -m skills.04_morphogenesis.vectorfield \
  --stage1 ./data/stage1.h5ad --stage2 ./data/stage2.h5ad
```

All skills accept `--help` for full parameter lists.

## Pipeline

```
Stage 0 (IO) ──→ Stage 1 (Alignment) ─→ Stage 2 (Reconstruction)
                                                    ↓
                                     interpolation    morphology
                                                    ↓
                                        Stage 3 (Morphogenesis)
                                        vectorfield ──→ feature
```

| Stage | Skills | Output |
|-------|--------|--------|
| **0: IO** | `manual_read`, `auto_read` | AnnData with `.obsm["spatial"]` |
| **1: Alignment** | `two_slice`, `multi_slice`, `two_stage_3d` | AnnData with `.obsm["align_spatial"]` |
| **2: Reconstruction** | `reconstruction`, `interpolation`, `morphology` | VTK models + interpolated expression |
| **3: Morphogenesis** | `vectorfield`, `feature` | Cell mappings, vector fields, GLM DEG tables |

## Supported Platforms (Stage 0)

| Platform | Reader | Input Format |
|----------|--------|-------------|
| MERFISH | `read_merfish` | `cell_by_gene.csv`, `cell_metadata.csv` |
| seqFISH | `read_seqfish` | CxG CSV, coordinates CSV, DAPI TIFF |
| Slide-seq | `read_slideseq` | Standard output directory |
| STARmap+ | `read_starmap_plus` | Expression/spatial CSVs |
| Stereo-seq | `read_stereoseq_bgi` / `read_stereoseq_agg` | GEM file |
| Xenium | `read_xenium` | 10x outs directory |
| Visium | `read_visium` | 10x Visium output |
| Visium HD | `read_visium_hd` | 10x Visium HD binned output |
| Open-ST | `read_openst` | Open-ST output |

## Project Structure

```
Spateo-Skills/
├── docs/
│   ├── README.md              # Project guide
│   ├── architecture/
│   │   ├── PYTHON-PATTERNS.md # Coding standards
│   │   └── CODEMAP.md         # File index
│   └── guides/                # How-to guides (future)
├── scripts/
│   └── sync-to-codex-plugin.sh # Plugin sync utility
├── skills/
│   ├── 00_shared/             # Shared infrastructure
│   ├── 00_env_setup/          # Environment setup
│   ├── 01_data_io/            # Stage 0: Data readers
│   ├── 02_slice_alignment/    # Stage 1: Alignment
│   ├── 03_3d_reconstruction/  # Stage 2: 3D models
│   └── 04_morphogenesis/      # Stage 3: Development
├── tests/                     # Test suite
│   ├── data_io/
│   ├── slice_alignment/
│   ├── reconstruction/
│   └── morphogenesis/
├── README.md                  # ← you are here
├── LICENSE
├── pyproject.toml
└── CHANGELOG.md
```

## Design Principles

- **Lazy imports** — importing a skill module does not load Spateo or execute
- **Parameterize everything** — no hardcoded paths, keys, or magic numbers
- **Validate at boundaries** — check file existence, array shapes, NaN/Inf
- **Composable** — each stage output feeds the next
- **CLI-ready** — every skill runs standalone from the command line

## Development

See `docs/README.md` for the full project guide, coding standards, and
instructions for adding new skills.
