#!/usr/bin/env python3
"""manual_read — Platform-specific spatial transcriptomics readers.

Wraps actual spateo IO functions: st.io.read_10x, st.io.read_bgi,
st.io.read_slideseq, st.io.read_nanostring, st.io.read_image.
Converts raw data into AnnData with .obsm["spatial"] populated.

Usage:
    python -m skills.01_data_io.manual_read --platform xenium --data-path ./data/xenium_outs/

Supported: MERFISH, seqFISH, Slide-seq, STARmap+, Stereo-seq, Xenium,
Visium, Visium HD, Open-ST.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Union

import numpy as np

LOGGER = logging.getLogger("data_io.manual_read")

PathLike = Union[str, Path]


def _as_str(path: PathLike) -> str:
    return str(Path(path))


def _import_spateo() -> Any:
    try:
        import spateo as st
        return st
    except ImportError as exc:
        raise ImportError(
            "spateo is required. Install via Spateo-Skills/00_env_setup/."
        ) from exc


def read_10x_visium(
    data_dir: PathLike,
    *,
    sample: Optional[str] = None,
    **kwargs: Any,
) -> Any:
    """Read 10x Visium output directory."""
    st = _import_spateo()
    return st.io.read_10x(_as_str(data_dir), sample=sample, **kwargs)


def read_10x_visium_hd(
    data_dir: PathLike,
    *,
    sample: Optional[str] = None,
    **kwargs: Any,
) -> Any:
    """Read 10x Visium HD binned output."""
    st = _import_spateo()
    return st.io.read_10x(_as_str(data_dir), sample=sample, **kwargs)


def read_xenium(
    data_dir: PathLike,
    **kwargs: Any,
) -> Any:
    """Read 10x Xenium output directory."""
    st = _import_spateo()
    return st.io.read_10x(_as_str(data_dir), **kwargs)


def read_stereoseq(
    gem_file: PathLike,
    *,
    binsize: int = 50,
    **kwargs: Any,
) -> Any:
    """Read Stereo-seq GEM file."""
    st = _import_spateo()
    return st.io.read_bgi(_as_str(gem_file), binsize=binsize, **kwargs)


def read_stereoseq_agg(
    gem_file: PathLike,
    image_file: PathLike,
    **kwargs: Any,
) -> Any:
    """Read Stereo-seq GEM + image (aggregated)."""
    st = _import_spateo()
    return st.io.read_bgi_agg(_as_str(gem_file), _as_str(image_file), **kwargs)


def read_slideseq(
    data_dir: PathLike,
    *,
    beads_path: Optional[str] = None,
    **kwargs: Any,
) -> Any:
    """Read Slide-seq output directory."""
    st = _import_spateo()
    if beads_path is None:
        beads_path = _as_str(Path(data_dir) / "BeadLocationsForR.csv")
    return st.io.read_slideseq(_as_str(data_dir), beads_path=beads_path, **kwargs)


def read_merfish(
    data_dir: PathLike,
    *,
    counts_file: str = "cell_by_gene.csv",
    meta_file: str = "cell_metadata.csv",
    **kwargs: Any,
) -> Any:
    """Read MERFISH from CSV files (gene matrix + cell metadata with coords)."""
    import anndata as ad
    import pandas as pd

    counts = pd.read_csv(Path(data_dir) / counts_file, index_col=0)
    meta = pd.read_csv(Path(data_dir) / meta_file, index_col=0)

    spatial_col_pairs = [
        ("X_center", "Y_center"), ("center_x", "center_y"),
        ("x", "y"), ("X", "Y"),
    ]
    for col_x, col_y in spatial_col_pairs:
        if col_x in meta.columns and col_y in meta.columns:
            spatial = meta[[col_x, col_y]].values
            break
    else:
        raise ValueError(f"No spatial columns found in {meta_file}")

    adata = ad.AnnData(X=counts.values, obs=meta, var=pd.DataFrame(index=counts.columns))
    adata.obsm["spatial"] = spatial
    return adata


def read_seqfish(
    data_dir: PathLike,
    *,
    counts_file: str,
    meta_file: str,
    load_images: bool = True,
    **kwargs: Any,
) -> Any:
    """Read seqFISH from CxG CSV + coordinates CSV."""
    import anndata as ad
    import pandas as pd

    counts = pd.read_csv(Path(data_dir) / counts_file, index_col=0)
    coords = pd.read_csv(Path(data_dir) / meta_file, index_col=0)

    x_col = next((c for c in coords.columns if c.lower().startswith("x")), None)
    y_col = next((c for c in coords.columns if c.lower().startswith("y")), None)
    if x_col and y_col:
        spatial = coords[[x_col, y_col]].values
    else:
        spatial = coords.iloc[:, :2].values

    adata = ad.AnnData(X=counts.values, obs=coords, var=pd.DataFrame(index=counts.columns))
    adata.obsm["spatial"] = spatial
    return adata


def read_starmap_plus(
    data_dir: PathLike,
    *,
    counts_file: str = "sagittal1processed_expression_pd.csv",
    meta_file: Optional[str] = None,
    spatial_file: Optional[str] = None,
    **kwargs: Any,
) -> Any:
    """Read STARmap+ from expression CSV + optional spatial metadata.

    Expression file: rows=spots, cols=genes.
    Spatial coords come from the spatial_file if provided, otherwise
    defaults are used.
    """
    import anndata as ad
    import pandas as pd

    counts = pd.read_csv(Path(data_dir) / counts_file, index_col=0)
    var_names = counts.columns.tolist()

    obs = {}
    if meta_file:
        meta = pd.read_csv(Path(data_dir) / meta_file, index_col=0)
        obs = meta

    spatial = None
    if spatial_file:
        spatial_df = pd.read_csv(Path(data_dir) / spatial_file)
        x_col = next((c for c in spatial_df.columns if c.lower() in ("x", "name")), None)
        y_col = next((c for c in spatial_df.columns if c.lower() == "y"), None)
        if x_col and y_col:
            spatial = spatial_df[[x_col, y_col]].values

    adata = ad.AnnData(X=counts.values, obs=pd.DataFrame(obs) if obs else None, var=pd.DataFrame(index=var_names))
    if spatial is not None:
        adata.obsm["spatial"] = spatial
    return adata


def read_openst(
    data_dir: PathLike,
    *,
    h5ad_file: Optional[str] = None,
    **kwargs: Any,
) -> Any:
    """Read Open-ST from h5ad or directory."""
    st = _import_spateo()
    if h5ad_file:
        return st.read(_as_str(Path(data_dir) / h5ad_file), **kwargs)
    return st.read(_as_str(data_dir), **kwargs)


READERS: Dict[str, Callable[..., Any]] = {
    "merfish": read_merfish,
    "seqfish": read_seqfish,
    "slideseq": read_slideseq,
    "starmap_plus": read_starmap_plus,
    "starmap-plus": read_starmap_plus,
    "stereoseq": read_stereoseq,
    "stereo-seq": read_stereoseq,
    "stereoseq_agg": read_stereoseq_agg,
    "stereo-seq-agg": read_stereoseq_agg,
    "xenium": read_xenium,
    "visium": read_10x_visium,
    "visium_hd": read_10x_visium_hd,
    "visium-hd": read_10x_visium_hd,
    "openst": read_openst,
    "open-st": read_openst,
}


def read_by_platform(platform: str, data_path: PathLike, **kwargs: Any) -> Any:
    """Dispatch to a platform-specific reader by name."""
    key = platform.lower().replace(" ", "_")
    if key not in READERS:
        raise ValueError(
            f"Unsupported platform: {platform!r}. Supported: {sorted(READERS)}"
        )
    return READERS[key](data_path, **kwargs)


__all__ = [
    "read_merfish",
    "read_seqfish",
    "read_slideseq",
    "read_starmap_plus",
    "read_stereoseq",
    "read_stereoseq_agg",
    "read_xenium",
    "read_10x_visium",
    "read_10x_visium_hd",
    "read_openst",
    "read_by_platform",
    "READERS",
]
