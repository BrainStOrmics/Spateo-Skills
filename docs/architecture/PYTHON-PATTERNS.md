# PYTHON-PATTERNS.md — Spateo Skill Module Standards

Read when writing or reviewing Spateo skill modules.

## Module Template

```python
#!/usr/bin/env python3
"""<skill_name> — Skill: <one-line description>.
Source notebook: <notebook_name>.ipynb
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

LOGGER = logging.getLogger("<logger_name>")

# Config & Result dataclasses
# Lazy imports (_import_spateo)
# Atomic functions (one responsibility)
# run_<name>_pipeline(config) -> result
# CLI (build_arg_parser, main)
```

## Rules

- `from __future__ import annotations` always
- Type hints on all signatures
- `**kwargs` passthrough for Spateo API flexibility
- `__all__` at bottom for public API
- No module-level side effects
- Config: `Path` for paths, `Optional[T]` for optional, defaults for safe runs
- Result: paths as `str` (JSON serializable), `notes: List[str]` for warnings
- Logging: `LOGGER` at module, `configure_logging(verbose)` in `main()`
- JSON output: `json.dumps(asdict(result), indent=2, ensure_ascii=False)`

## Naming

- Functions: `snake_case`, verb-first (`construct_point_cloud`)
- Config: `<SkillName>Config`
- Result: `<SkillName>Result`
- AnnData keys: `"align_spatial"`, `"3d_align_spatial"`, `"VecFld_morpho"`

## Validation (at boundaries only)

- File/directory existence
- Spatial key in `adata.obsm`
- Array shapes (`coords.shape[1] >= 3` for 3D)
- No NaN/Inf in critical arrays
- Gene names in `adata.var_names`

## Error Handling

- Let Spateo raise; wrap with context if helpful
- Pipeline catches in `main()`, returns exit code 1
- Multi-method skills: catch per-method, store error in result, continue

## Shared Infrastructure

All skills import from `00_shared/`:

- `from ..shared.lazy_imports import _import_spateo, _import_spateo_plus`
- `from ..shared.data_helpers import (configure_logging, infer_device, load_h5ad, validate_spatial, LOGGER, PathLike)`

Relative imports work because all skills are under `skills/` as a package.
