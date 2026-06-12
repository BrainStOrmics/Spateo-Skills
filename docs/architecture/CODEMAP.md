# CODEMAP.md — Spateo-Skills File Index

## Entry Points
| File | Purpose |
|------|---------|
| `docs/README.md` | Project overview, quickstart, pipeline docs |
| `docs/architecture/PYTHON-PATTERNS.md` | Coding standards for all skill modules |

## Infrastructure
| File | Purpose |
|------|---------|
| `skills/00_shared/__init__.py` | Package init |
| `skills/00_shared/lazy_imports.py` | `_import_spateo()`, `_import_spateo_plus()`, `_import_pyvista()` |
| `skills/00_shared/data_helpers.py` | `load_h5ad()`, `validate_spatial()`, `infer_device()`, `configure_logging()` |
| `skills/00_env_setup/setup.sh` | Interactive env installer — auto-detects conda/uv, prompts choice |
| `skills/00_env_setup/environment.yml` | conda-forge spec (pip fallback for spateo-release) |
| `skills/00_env_setup/uv_requirements.txt` | Flat requirements for `uv pip install` |
| `skills/00_env_setup/uv_pyproject.toml` | pyproject.toml with `[tool.uv.sources]` git deps |
| `skills/00_env_setup/skill.md` | conda + uv setup docs with troubleshooting |

## Pipeline Orchestrator
| File | Purpose |
|------|---------|
| `skills/00_pipeline/__init__.py` | Package init |
| `skills/00_pipeline/run_all.py` | Full pipeline: IO → Alignment → Reconstruction → Morphogenesis |
| `skills/00_pipeline/SKILL.md` | Pipeline docs with per-stage skip flags |

## Stage 0: Data IO
| File | Purpose |
|------|---------|
| `skills/01_data_io/manual_read.py` | 10 platform-specific readers + `read_by_platform()` dispatcher |
| `skills/01_data_io/auto_read.py` | Auto-detect + batch readers |
| `skills/01_data_io/SKILL.md` | Platform readers reference |

## Stage 1: Slice Alignment
| File | Purpose |
|------|---------|
| `skills/02_slice_alignment/two_slice.py` | Two-slice `morpho_align` pipeline |
| `skills/02_slice_alignment/multi_slice.py` | Continuous multi-slice alignment + 3D concat |
| `skills/02_slice_alignment/two_stage_3d.py` | Two-stage 3D model `morpho_align_ref` |
| `skills/02_slice_alignment/SKILL.md` | Alignment reference (morpho_align / morpho_align_ref) |

## Stage 2: 3D Reconstruction
| File | Purpose |
|------|---------|
| `skills/03_3d_reconstruction/reconstruction.py` | Point cloud, surface, cell mesh, voxel, subtype models |
| `skills/03_3d_reconstruction/interpolation.py` | VTK/GP/kernel/deep gene expression interpolation |
| `skills/03_3d_reconstruction/morphology.py` | KDE density + morphology metrics |
| `skills/03_3d_reconstruction/SKILL.md` | 3D models and interpolation reference |

## Stage 3: Morphogenesis
| File | Purpose |
|------|---------|
| `skills/04_morphogenesis/vectorfield.py` | Cell mapping, `morphofield_sparsevfc`, `morphopath` |
| `skills/04_morphogenesis/feature.py` | velocity/acceleration/curvature/curl/torsion/jacobian + GLM DEG |
| `skills/04_morphogenesis/SKILL.md` | Vector field and GLM DEG reference |

## Tests
| Directory | Coverage |
|-----------|----------|
| `tests/data_io/` | Stage 0 readers |
| `tests/slice_alignment/` | Stage 1 alignment |
| `tests/reconstruction/` | Stage 2 reconstruction/interpolation/morphology |
| `tests/morphogenesis/` | Stage 3 vectorfield/feature |

## Scripts
| File | Purpose |
|------|---------|
| `scripts/sync-to-codex-plugin.sh` | Sync skills to external plugin registries |
