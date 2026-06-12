# Spateo-Skills — Project Guide

## What This Is
Standardized, parameterized Python skill modules for the Spateo spatial
transcriptomics pipeline. Each module converts a notebook workflow into a
reusable, CLI-runnable skill with typed Config/Result dataclasses, lazy imports,
and boundary validation.

## Architecture

```
skills/
├── 00_shared/          # Shared infrastructure (lazy imports, data helpers)
├── 00_env_setup/       # Environment setup (conda/venv, requirements)
├── 01_data_io/         # Stage 0: platform-specific + auto-detect readers
├── 02_slice_alignment/ # Stage 1: 2D→3D coordinate registration
├── 03_3d_reconstruction/ # Stage 2: point cloud, mesh, interpolation, morphology
└── 04_morphogenesis/   # Stage 3: vector field, trajectory, GLM DEG
```

Full pipeline: `IO → Alignment → Reconstruction → Morphogenesis`

## Coding Standards

Every skill module follows `docs/architecture/PYTHON-PATTERNS.md`:

1. `from __future__ import annotations`
2. Type hints on all signatures
3. `**kwargs` passthrough for Spateo API flexibility
4. `__all__` at bottom for public API
5. Lazy imports — no module-level side effects
6. `@dataclass <SkillName>Config` — all parameters with defaults
7. `@dataclass <SkillName>Result` — paths as `str`, `notes: List[str]`
8. Atomic functions (one responsibility, typed)
9. `run_<name>_pipeline(config) -> result` — end-to-end orchestrator
10. Optional `build_arg_parser()`, `main()` for CLI
11. Validation at boundaries only

## Quickstart

### Environment

```bash
bash skills/00_env_setup/setup.sh
# or manually:
conda create -n spateo python=3.10 -y && conda activate spateo
pip install -r skills/00_env_setup/requirements.txt
```

### Run a Skill

```bash
# Stage 0: Read data
python -m skills.01_data_io.manual_read --platform xenium --data-path /path/to/data

# Stage 1: Align two slices
python -m skills.02_slice_alignment.two_slice --slice1 a.h5ad --slice2 b.h5ad

# Stage 2: 3D reconstruction
python -m skills.03_3d_reconstruction.reconstruction --input aligned.h5ad --out-dir ./output

# Stage 3: Morphogenesis
python -m skills.04_morphogenesis.vectorfield --stage1 s1.h5ad --stage2 s2.h5ad
```

### Test Data

Test data lives in the parent project at `../data/skills_data/`. Skills reference
it via relative paths in their default Config values.

## Development

### Adding a New Skill

1. Create directory: `skills/NN_name/`
2. Add `__init__.py`, `<skill>.py`
3. Follow PYTHON-PATTERNS.md template
4. Add tests in `tests/<category>/`
5. Update docs/architecture/CODEMAP.md
6. Update CHANGELOG.md

### Running Tests

```bash
pytest tests/ -v
# or individual category:
pytest tests/data_io/ -v
```

## Pipeline Dependencies

```
Stage 0 (IO) ──→ Stage 1 (Alignment) ─→ Stage 2 (Reconstruction)
                                                    ↓
                                 interpolation   morphology
                                                    ↓
                                        Stage 3 (Morphogenesis)
                                        vectorfield → feature
```

Each stage consumes `.h5ad` from the previous and produces annotated `.h5ad` + `.vtk` models.
