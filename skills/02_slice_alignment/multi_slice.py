"""Continuous multi-slice spatial alignment and 3D reconstruction skill.

Splits data by z-coordinate, aligns consecutive pairs via
morpho_align_transformation, and concatenates into 3D.

Usage:
    python -m skills.02_slice_alignment.multi_slice \
        --input ./data/slices/ --outdir ./output/

Source: drosophila_Alignment(1).ipynb
Data:   ../../data/skills_data/1.slices alignment/multi_slices/
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import numpy as np

from skills.shared.data_helpers import (
    LOGGER,
    PathLike,
    configure_logging,
    infer_device,
    load_h5ad,
    validate_spatial,
)

_DATA_ROOT = Path(__file__).resolve().parent.parent.parent / "data" / "skills_data" / "1.slices alignment" / "multi_slices"


# ---------------------------------------------------------------------------
# Config / Result
# ---------------------------------------------------------------------------


@dataclass
class MultiSliceConfig:
    adata_path: str = str(_DATA_ROOT / "E16-18h_a_count_normal_stereoseq.h5ad")
    preprocess: bool = True
    recipe: str = "pearson_residuals"
    spatial_key: str = "spatial"
    counts_layer: str = "counts"
    n_top_genes: int = 3000
    run_pca: bool = True
    xyz_obs_keys: tuple[str, str, str] = ("new_x", "new_y", "new_z")
    label_key: str = "annotation"
    spatial_3d_key: str = "spatial_3D"
    spatial_2d_key: str = "spatial_2D"
    key_added: str = "align_spatial"
    output_3d_key: str = "aligned_spatial_3D"
    initialize_rigid: bool = False
    rep_layer: str = "X_pca"
    rep_field: str = "obsm"
    dissimilarity: str = "cos"
    device: Optional[str] = None
    cuda_visible_devices: Optional[str] = None
    outdir: Optional[str] = None
    transform_kwargs: dict = field(default_factory=dict)


@dataclass
class MultiSliceResult:
    adata_path: str = ""
    n_slices: int = 0
    n_obs_total: int = 0
    aligned_2d_key: str = ""
    aligned_3d_key: str = ""
    device: str = ""
    outdir: str = ""


# ---------------------------------------------------------------------------
# Atomic helpers
# ---------------------------------------------------------------------------


def load_adata(path: PathLike) -> Any:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"adata not found: {p}")
    return load_h5ad(p)


def preprocess_adata(
    adata: Any,
    *,
    recipe: str = "pearson_residuals",
    spatial_key: str = "spatial",
    counts_layer: str = "counts",
    n_top_genes: int = 3000,
    run_pca: bool = True,
    **kwargs: Any,
) -> Any:
    import spateo as st

    preprocess_spatial(
        adata,
        recipe=recipe,
        spatial_key=spatial_key,
        counts_layer=counts_layer,
        n_top_genes=n_top_genes,
        run_pca=run_pca,
        **kwargs,
    )
    return adata


def split_by_z(
    adata: Any,
    *,
    xyz_obs_keys: tuple[str, str, str] = ("new_x", "new_y", "new_z"),
    spatial_3d_key: str = "spatial_3D",
    spatial_2d_key: str = "spatial_2D",
    sort_z: bool = True,
    initialize_rigid: bool = False,
) -> Sequence[Any]:
    adata.obsm[spatial_3d_key] = adata.obs[list(xyz_obs_keys)].values.copy()
    z_height = np.unique(adata.obsm[spatial_3d_key][:, 2])
    if sort_z:
        z_height = np.sort(z_height)
    slices = [adata[adata.obsm[spatial_3d_key][:, 2] == z].copy() for z in z_height]
    for i, slc in enumerate(slices):
        slc.obsm[spatial_2d_key] = slc.obsm[spatial_3d_key][:, :2].copy()
        if initialize_rigid and i > 0:
            import spateo as st

            st.align.rigid_transformation(slc, spatial_key=spatial_2d_key, key_added=spatial_2d_key)
    return slices


def compute_multi_slice_transformation(
    slices: Sequence[Any],
    *,
    spatial_key: str = "spatial_2D",
    key_added: str = "align_spatial",
    device: Optional[str] = None,
    rep_layer: str = "X_pca",
    rep_field: str = "obsm",
    dissimilarity: str = "cos",
    verbose: bool = False,
    **kwargs: Any,
) -> Any:
    import spateo as st

    return st.align.morpho_align_transformation(
        models=list(slices),
        spatial_key=spatial_key,
        key_added=key_added,
        device=device or infer_device(),
        verbose=verbose,
        rep_layer=rep_layer,
        rep_field=rep_field,
        dissimilarity=dissimilarity,
        **kwargs,
    )


def apply_multi_slice_transformation(
    slices: Sequence[Any],
    transformation: Any,
    *,
    spatial_key: str = "spatial_2D",
    key_added: str = "align_spatial",
    **kwargs: Any,
) -> Sequence[Any]:
    import spateo as st

    return st.align.morpho_align_apply_transformation(
        models=list(slices),
        spatial_key=spatial_key,
        key_added=key_added,
        transformation=transformation,
        **kwargs,
    )


def concatenate_aligned_slices_to_3d(
    aligned_slices: Sequence[Any],
    *,
    aligned_2d_key: str = "align_spatial",
    original_3d_key: str = "spatial_3D",
    output_3d_key: str = "aligned_spatial_3D",
) -> Any:
    import anndata as ad

    aligned_adata = ad.concat(list(aligned_slices))
    z = np.asarray(aligned_adata.obsm[original_3d_key][:, 2])[:, None]
    aligned_adata.obsm[output_3d_key] = np.concatenate(
        [aligned_adata.obsm[aligned_2d_key], z], axis=1
    )
    return aligned_adata


def construct_aligned_point_cloud(
    aligned_adata: Any,
    *,
    spatial_key: str = "aligned_spatial_3D",
    groupby: str = "annotation",
    key_added: str = "tissue",
    colormap: Any = "rainbow",
    **kwargs: Any,
) -> tuple[Any, Any]:
    import spateo as st

    return st.tdr.construct_pc(
        adata=aligned_adata,
        spatial_key=spatial_key,
        groupby=groupby,
        key_added=key_added,
        colormap=colormap,
        **kwargs,
    )


def plot_slices(
    slices: Sequence[Any],
    *,
    label_key: str = "annotation",
    spatial_key: str = "spatial_2D",
    height: int = 2,
    ncols: int = 5,
    **kwargs: Any,
) -> Any:
    import spateo as st

    return st.pl.slices_2d(
        slices=list(slices),
        label_key=label_key,
        spatial_key=spatial_key,
        height=height,
        center_coordinate=False,
        show_legend=True,
        ncols=ncols,
        **kwargs,
    )


def plot_aligned_point_cloud(
    aligned_pc: Any,
    *,
    key: str = "tissue",
    **kwargs: Any,
) -> Any:
    import spateo as st

    defaults = dict(model_style="points", model_size=8, show_axes=True, jupyter="static")
    defaults.update(kwargs)
    return st.pl.three_d_plot(model=aligned_pc, key=key, **defaults)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def run_multi_slice_pipeline(config: MultiSliceConfig) -> MultiSliceResult:
    device = infer_device(config.cuda_visible_devices) if config.device is None else config.device
    LOGGER.info("Multi-slice alignment: device=%s", device)

    adata = load_adata(config.adata_path)
    LOGGER.info("Loaded adata: n_obs=%d", adata.n_obs)

    if config.preprocess:
        preprocess_adata(
            adata,
            recipe=config.recipe,
            spatial_key=config.spatial_key,
            counts_layer=config.counts_layer,
            n_top_genes=config.n_top_genes,
            run_pca=config.run_pca,
        )

    slices = split_by_z(
        adata,
        xyz_obs_keys=config.xyz_obs_keys,
        spatial_3d_key=config.spatial_3d_key,
        spatial_2d_key=config.spatial_2d_key,
        initialize_rigid=config.initialize_rigid,
    )
    LOGGER.info("Split into %d z-slices", len(slices))
    for i, slc in enumerate(slices):
        validate_spatial(slc, config.spatial_2d_key)

    transformation = compute_multi_slice_transformation(
        slices,
        spatial_key=config.spatial_2d_key,
        key_added=config.key_added,
        device=device,
        rep_layer=config.rep_layer,
        rep_field=config.rep_field,
        dissimilarity=config.dissimilarity,
        **dict(config.transform_kwargs),
    )
    aligned_slices = apply_multi_slice_transformation(
        slices,
        transformation,
        spatial_key=config.spatial_2d_key,
        key_added=config.key_added,
    )
    aligned_adata = concatenate_aligned_slices_to_3d(
        aligned_slices,
        aligned_2d_key=config.key_added,
        original_3d_key=config.spatial_3d_key,
        output_3d_key=config.output_3d_key,
    )
    aligned_pc, _cmap = construct_aligned_point_cloud(
        aligned_adata,
        spatial_key=config.output_3d_key,
        groupby=config.label_key,
    )

    outdir = ""
    if config.outdir:
        out = Path(config.outdir)
        out.mkdir(parents=True, exist_ok=True)
        outdir = str(out)

    return MultiSliceResult(
        adata_path=str(config.adata_path),
        n_slices=len(slices),
        n_obs_total=adata.n_obs,
        aligned_2d_key=config.key_added,
        aligned_3d_key=config.output_3d_key,
        device=device,
        outdir=outdir,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Continuous multi-slice spatial alignment")
    p.add_argument("--adata", type=str, default=None)
    p.add_argument("--no-preprocess", action="store_true")
    p.add_argument("--recipe", type=str, default="pearson_residuals")
    p.add_argument("--spatial-key", type=str, default="spatial")
    p.add_argument("--spatial-2d-key", type=str, default="spatial_2D")
    p.add_argument("--spatial-3d-key", type=str, default="spatial_3D")
    p.add_argument("--key-added", type=str, default="align_spatial")
    p.add_argument("--output-3d-key", type=str, default="aligned_spatial_3D")
    p.add_argument("--label-key", type=str, default="annotation")
    p.add_argument("--initialize-rigid", action="store_true")
    p.add_argument("--n-top-genes", type=int, default=3000)
    p.add_argument("--device", type=str, default=None)
    p.add_argument("--outdir", type=str, default=None)
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_arg_parser().parse_args(argv)
    configure_logging(args.verbose)

    cfg = MultiSliceConfig(
        preprocess=not args.no_preprocess,
        recipe=args.recipe,
        spatial_key=args.spatial_key,
        spatial_2d_key=args.spatial_2d_key,
        spatial_3d_key=args.spatial_3d_key,
        key_added=args.key_added,
        output_3d_key=args.output_3d_key,
        label_key=args.label_key,
        initialize_rigid=args.initialize_rigid,
        n_top_genes=args.n_top_genes,
        device=args.device,
        outdir=args.outdir,
    )
    if args.adata:
        cfg.adata_path = args.adata

    result = run_multi_slice_pipeline(cfg)
    print(json.dumps(asdict(result), indent=2, ensure_ascii=False))


__all__ = [
    "MultiSliceConfig",
    "MultiSliceResult",
    "load_adata",
    "preprocess_adata",
    "split_by_z",
    "compute_multi_slice_transformation",
    "apply_multi_slice_transformation",
    "concatenate_aligned_slices_to_3d",
    "construct_aligned_point_cloud",
    "plot_slices",
    "plot_aligned_point_cloud",
    "run_multi_slice_pipeline",
    "build_arg_parser",
    "main",
]
