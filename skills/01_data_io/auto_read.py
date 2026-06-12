#!/usr/bin/env python3
"""auto_read — Auto-detect platform and read spatial transcriptomics data.

Detects the spatial platform (Xenium, Visium, MERFISH, etc.) from the
data directory structure, then reads data into AnnData with
``.obsm["spatial"]`` populated.

Usage:
    python -m skills.01_data_io.auto_read --data-path ./data/unknown/

Source notebook: autoread_skills.ipynb
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Iterable, Union

LOGGER = logging.getLogger("data_io.auto_read")

PathLike = Union[str, Path]


def _import_spateo() -> Any:
    try:
        import spateo as st  # type: ignore
        return st
    except ImportError as exc:
        raise ImportError(
            "spateo is required. Install via Spateo-Skills/00_env_setup/."
        ) from exc


def _load_auto_io() -> tuple[Any, Any]:
    try:
        from spateo.io.protocol_io.spatial.auto import (  # type: ignore
            detect_spatial_technology,
            read_auto_spatial,
        )
        return detect_spatial_technology, read_auto_spatial
    except Exception as exc:
        raise ImportError(
            "spateo.io.protocol_io.spatial.auto is required for auto-detection."
        ) from exc


def detect_spatial_platform(data_dir: PathLike, **detect_kwargs: Any) -> Any:
    """Auto-detect the spatial transcriptomics platform for a data directory."""
    detect_spatial_technology, _ = _load_auto_io()
    return detect_spatial_technology(str(Path(data_dir)), **detect_kwargs)


def read_auto_spatial_data(
    data_dir: PathLike,
    *,
    return_match: bool = False,
    **read_kwargs: Any,
) -> Any:
    """Auto-detect platform and read spatial data into AnnData.

    Returns ``(adata, match)`` when *return_match* is True, else ``adata``.
    """
    detect_spatial_technology, read_auto_spatial = _load_auto_io()
    data_dir_str = str(Path(data_dir))
    match = detect_spatial_technology(data_dir_str)
    adata = read_auto_spatial(data_dir_str, **read_kwargs)
    return (adata, match) if return_match else adata


def read_many_auto(
    data_dirs: Iterable[PathLike],
    *,
    return_matches: bool = True,
    stop_on_error: bool = False,
    **read_kwargs: Any,
) -> Dict[str, Any]:
    """Batch auto-read multiple spatial data directories.

    On failure per directory: stores the exception in the result dict
    unless *stop_on_error* is True.
    """
    results: Dict[str, Any] = {}
    for data_dir in data_dirs:
        key = str(Path(data_dir))
        try:
            if return_matches:
                adata, match = read_auto_spatial_data(
                    key, return_match=True, **read_kwargs
                )
                results[key] = {"adata": adata, "match": match, "error": None}
            else:
                results[key] = read_auto_spatial_data(
                    key, return_match=False, **read_kwargs
                )
        except Exception as exc:
            if stop_on_error:
                raise
            LOGGER.warning("Failed to read %s: %s", key, exc)
            results[key] = {"adata": None, "match": None, "error": exc}
    return results


__all__ = [
    "detect_spatial_platform",
    "read_auto_spatial_data",
    "read_many_auto",
]
