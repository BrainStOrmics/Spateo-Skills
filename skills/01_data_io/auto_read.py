#!/usr/bin/env python3
"""auto_read — Auto-detect platform and read spatial transcriptomics data.

Detects the spatial platform by inspecting directory contents (file names,
structure), then dispatches to the appropriate reader from manual_read.

Usage:
    python -m skills.01_data_io.auto_read --data-path ./data/unknown_dataset/
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Iterable, Union

LOGGER = logging.getLogger("data_io.auto_read")

PathLike = Union[str, Path]

PLATFORM_SIGNATURES = {
    "xenium": [
        "cell_boundaries.csv", "cell_boundaries.csv.gz",
        "cell_boundaries.parquet", "cells.csv", "cells.csv.gz",
        "cells.parquet", "gene_panel.json",
    ],
    "visium": [
        "spatial/scalefactors_json.json",
        "spatial/tissue_positions_list.csv",
    ],
    "visium_hd": [
        "spatial/tissue_positions.parquet",
        "binned_outputs",
    ],
    "merfish": [
        "cell_by_gene.csv", "cell_metadata.csv",
    ],
    "seqfish": [
        "CxG", "_CxG_", "section",
    ],
    "stereoseq": [
        ".gem", ".gem.txt", ".bin",
    ],
    "slideseq": [
        "MappedDGEForR.csv", "BeadLocationsForR.csv",
    ],
    "starmap_plus": [
        "processed_expression_pd.csv", "_spatial.csv",
    ],
    "openst": [
        ".h5ad",
    ],
}


def _files_in_dir(path: Path) -> set[str]:
    """Return set of file and subdir names in *path*."""
    return {p.name for p in path.iterdir()} if path.is_dir() else set()


def _file_exists(path: Path, subpath: str) -> bool:
    return (path / subpath).exists()


def detect_spatial_platform(data_dir: PathLike, **kwargs: Any) -> str:
    """Auto-detect the spatial transcriptomics platform for a data directory.

    Checks file names and sub-directory patterns against known signatures.
    Returns the platform key (e.g. "xenium", "visium", "merfish").
    """
    p = Path(data_dir)
    files = _files_in_dir(p)
    file_list = " ".join(files).lower()

    for platform, signatures in PLATFORM_SIGNATURES.items():
        for sig in signatures:
            sig_lower = sig.lower()
            if "/" in sig:
                if _file_exists(p, sig):
                    return platform
            else:
                if any(sig_lower in f.lower() for f in files):
                    return platform

    raise ValueError(
        f"Cannot detect platform for {data_dir!r}. "
        f"Recognized signatures: {list(PLATFORM_SIGNATURES.keys())}"
    )


def read_auto_spatial_data(
    data_dir: PathLike,
    *,
    return_match: bool = False,
    **read_kwargs: Any,
) -> Any:
    """Auto-detect platform and read spatial data into AnnData.

    Returns (adata, platform) when return_match is True, else adata.
    """
    from .manual_read import read_by_platform

    platform = detect_spatial_platform(data_dir)
    LOGGER.info("Auto-detected platform: %s", platform)
    adata = read_by_platform(platform, data_dir, **read_kwargs)
    return (adata, platform) if return_match else adata


def read_many_auto(
    data_dirs: Iterable[PathLike],
    *,
    return_matches: bool = True,
    stop_on_error: bool = False,
    **read_kwargs: Any,
) -> Dict[str, Any]:
    """Batch auto-read multiple spatial data directories.

    On failure per directory: stores the exception in the result dict
    unless stop_on_error is True.
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
