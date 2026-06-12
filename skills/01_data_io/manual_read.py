#!/usr/bin/env python3
"""manual_read — Skill: platform-specific spatial transcriptomics readers.
Source notebooks: merfish(1), seqfish(1), slideseq(1), starmap_plus(1),
stereoseq(1), xenium(1), visium(1), visium_hd(2), openst.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Union

LOGGER = logging.getLogger("data_io.manual_read")

PathLike = Union[str, Path]


def _as_str(path: PathLike) -> str:
    return str(Path(path))


def _import_spateo() -> Any:
    try:
        import spateo as st  # type: ignore
        return st
    except ImportError as exc:
        raise ImportError(
            "spateo is required. Install via Spateo-Skills/00_env_setup/."
        ) from exc


def _spateo_protocol_io() -> Any:
    try:
        import spateo.io.protocol_io as pio  # type: ignore
        return pio
    except Exception as exc:
        raise ImportError(
            "spateo.io.protocol_io is required for this reader."
        ) from exc


def _spateo_plus_io_or_spateo_io() -> Any:
    try:
        import spateo_plus as sp  # type: ignore
        return sp.io
    except Exception:
        try:
            import spateo as st  # type: ignore
            return st.io
        except Exception as exc:
            raise ImportError(
                "spateo_plus or spateo with io reader is required."
            ) from exc


def read_merfish(
    data_dir: PathLike,
    *,
    counts_file: str = "cell_by_gene.csv",
    meta_file: str = "cell_metadata.csv",
    **kwargs: Any,
) -> Any:
    pio = _spateo_protocol_io()
    return pio.spatial.read_merfish(
        _as_str(data_dir), counts_file=counts_file, meta_file=meta_file, **kwargs
    )


def read_seqfish(
    data_dir: PathLike,
    *,
    counts_file: str,
    meta_file: str,
    load_images: bool = True,
    **kwargs: Any,
) -> Any:
    pio = _spateo_protocol_io()
    return pio.spatial.read_seqfish(
        _as_str(data_dir),
        counts_file=counts_file,
        meta_file=meta_file,
        load_images=load_images,
        **kwargs,
    )


def read_slideseq(
    data_dir: PathLike,
    *,
    load_images: bool = True,
    **kwargs: Any,
) -> Any:
    pio = _spateo_protocol_io()
    return pio.spatial.read_slideseq(
        _as_str(data_dir), load_images=load_images, **kwargs
    )


def read_starmap_plus(
    data_dir: PathLike,
    *,
    counts_file: str = "sagittal1processed_expression_pd.csv",
    meta_file: str = "sagittal1_spot_meta.csv",
    spatial_file: str = "sagittal1_spatial.csv",
    **kwargs: Any,
) -> Any:
    pio = _spateo_protocol_io()
    return pio.spatial.read_starmap_plus(
        _as_str(data_dir),
        counts_file=counts_file,
        meta_file=meta_file,
        spatial_file=spatial_file,
        **kwargs,
    )


def read_stereoseq(
    gem_file: PathLike,
    *,
    binsize: int = 50,
    **kwargs: Any,
) -> Any:
    pio = _spateo_protocol_io()
    return pio.spatial.read_bgi(
        _as_str(gem_file), binsize=binsize, **kwargs
    )


def read_stereoseq_agg(
    gem_file: PathLike,
    image_file: PathLike,
    **kwargs: Any,
) -> Any:
    pio = _spateo_protocol_io()
    return pio.spatial.read_bgi_agg(
        _as_str(gem_file), _as_str(image_file), **kwargs
    )


def read_xenium(data_dir: PathLike, **kwargs: Any) -> Any:
    io = _spateo_plus_io_or_spateo_io()
    return io.read_xenium(_as_str(data_dir), **kwargs)


def read_visium(data_dir: PathLike, **kwargs: Any) -> Any:
    io = _spateo_plus_io_or_spateo_io()
    return io.read_visium(_as_str(data_dir), **kwargs)


def read_visium_hd(data_dir: PathLike, **kwargs: Any) -> Any:
    io = _spateo_plus_io_or_spateo_io()
    return io.read_visium_hd(_as_str(data_dir), **kwargs)


def read_openst(
    data_dir: PathLike,
    *,
    h5ad_file: Optional[str] = None,
    **kwargs: Any,
) -> Any:
    io = _spateo_plus_io_or_spateo_io()
    if h5ad_file is None:
        return io.read_openst(_as_str(data_dir), **kwargs)
    return io.read_openst(_as_str(data_dir), h5ad_file=h5ad_file, **kwargs)


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
    "visium": read_visium,
    "visium_hd": read_visium_hd,
    "visium-hd": read_visium_hd,
    "openst": read_openst,
    "open-st": read_openst,
}


def read_by_platform(platform: str, data_path: PathLike, **kwargs: Any) -> Any:
    """Dispatch to a platform-specific reader by name.

    Examples
    --------
    >>> adata = read_by_platform("slideseq", "../../data/Slideseq", load_images=False)
    >>> adata = read_by_platform("xenium", "../../data/xenium_outs")
    """
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
    "read_visium",
    "read_visium_hd",
    "read_openst",
    "read_by_platform",
    "READERS",
]
