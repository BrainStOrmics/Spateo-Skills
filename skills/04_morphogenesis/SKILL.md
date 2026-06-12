---
name: spateo-morphogenesis
description: Vector field inference, trajectory analysis, GLM differential expression
---

# Stage 3: Morphogenesis

Infer cell state transitions and vector fields between developmental stages, compute trajectory paths, and identify differentially expressed genes.

## Quick Start

```bash
# Vector field + cell mapping between two stages
python -m skills.04_morphogenesis.vectorfield \
  --stage1 ./data/stage1.h5ad --stage2 ./data/stage2.h5ad

# Feature extraction (velocity, acceleration, curvature) + GLM DEG
python -m skills.04_morphogenesis.feature \
  --input ./data/vectorfield.h5ad --out-dir ./output
```

## Modules

### `vectorfield.py` — Cell mapping & vector field

Computes cell state transitions between two time points using sparse VFC (Vector Field Corrector).

```python
from skills.04_morphogenesis.vectorfield import (
    VectorFieldConfig, run_vectorfield_pipeline,
)

config = VectorFieldConfig(
    stage1_path="./data/stage1.h5ad",
    stage2_path="./data/stage2.h5ad",
    groupby="celltype",
    method="sparsevfc",
)
result = run_vectorfield_pipeline(config)
```

**Key APIs:** `st.tools.morphofield_sparsevfc`, `st.tools.morphopath`, cell mapping via optimal transport (POT).

### `feature.py` — Feature extraction + GLM DEG

Extracts dynamic features from vector fields and identifies differentially expressed genes.

```python
from skills.04_morphogenesis.feature import (
    FeatureConfig, run_feature_pipeline,
)

config = FeatureConfig(
    input_path="./data/vectorfield.h5ad",
    features=["velocity", "acceleration", "curvature", "curl", "torsion", "jacobian"],
    run_glm_deg=True,
    deg_groupby="celltype",
)
```

**Features:**

| Feature | Description |
|---------|-------------|
| `velocity` | Cell state change speed |
| `acceleration` | Rate of velocity change |
| `curvature` | Path bending |
| `curl` | Local rotation in vector field |
| `torsion` | 3D twisting of paths |
| `jacobian` | Local divergence/convergence |

## Input/Output

**Input (vectorfield):** Two AnnData objects from consecutive time stages.
**Input (feature):** AnnData with vector field computed.
**Output:** AnnData with vector field layers, GLM DEG tables (CSV), trajectory paths.

## Key Parameters

| Parameter | Purpose | Default |
|-----------|---------|---------|
| `stage1_path`, `stage2_path` | Input stage paths | — |
| `groupby` | Cell grouping column | `"celltype"` |
| `method` | Vector field method | `"sparsevfc"` |
| `features` | Features to compute | All |
| `run_glm_deg` | Run GLM differential expression | `True` |
| `deg_threshold` | Significance threshold | 0.05 |

## Troubleshooting

| Problem | Fix |
|---------|-----|
| POT import error | `conda install -c conda-forge pot` |
| Vector field diverges | Check stage1/stage2 have same gene set |
| GLM fails on sparse data | Filter low-expression genes first |
| CUDA OOM | Set `device="cpu"` in config |
