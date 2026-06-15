"""Common data helpers: load AnnData, validate spatial coords, device detection, logging."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np

LOGGER = logging.getLogger("spateo_shared")

PathLike = Union[str, Path]


def configure_logging(verbose: bool = False) -> None:
    """Configure root logger for all skill modules."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def infer_device(cuda_visible_devices: Optional[str] = None) -> str:
    """Return 'cuda' if GPU available, else 'cpu'."""
    if cuda_visible_devices is not None:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(cuda_visible_devices)
    try:
        import torch

        return "cuda" if torch.cuda.is_available() else "cpu"
    except Exception:
        return "cpu"


def load_h5ad(path: PathLike) -> Any:
    """Load AnnData via spateo.read_h5ad."""
    from .lazy_imports import _import_spateo

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"h5ad not found: {path}")
    st = _import_spateo()
    LOGGER.info("Reading AnnData: %s", path)
    return st.read_h5ad(str(path))


def validate_spatial(
    adata: Any,
    spatial_key: str,
    *,
    require_3d: bool = False,
) -> np.ndarray:
    """Validate spatial coordinates exist and are well-formed.

    Returns the coordinate array (n_obs x dim).
    """
    if spatial_key not in adata.obsm:
        raise KeyError(
            f"spatial key {spatial_key!r} not in adata.obsm. "
            f"Available: {list(adata.obsm.keys())}"
        )
    coords = np.asarray(adata.obsm[spatial_key])
    if coords.ndim != 2:
        raise ValueError(f"Spatial coords must be 2D, got shape {coords.shape}")
    if require_3d and coords.shape[1] < 3:
        raise ValueError(
            f"Need at least 3 columns for 3D spatial, got shape {coords.shape}"
        )
    if not np.isfinite(coords).all():
        raise ValueError(f"Spatial coords contain NaN/Inf at key {spatial_key!r}")
    return coords


__all__ = [
    "configure_logging",
    "infer_device",
    "load_h5ad",
    "validate_spatial",
    "LOGGER",
    "PathLike",
]
