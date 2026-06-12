"""Two-slice spatial alignment skill.

Source: cs13_Alignment(1).ipynb
Data:   ../../data/skills_data/1.slices alignment/2_slices/
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import numpy as np

from ..shared.data_helpers import (
    LOGGER,
    PathLike,
    configure_logging,
    infer_device,
    load_h5ad,
    validate_spatial,
)

_DATA_ROOT = Path(__file__).resolve().parent.parent.parent / "data" / "skills_data" / "1.slices alignment" / "2_slices"


# ---------------------------------------------------------------------------
# Config / Result
# ---------------------------------------------------------------------------


@dataclass
class TwoSliceConfig:
    slice1_path: str = str(_DATA_ROOT / "adata_109_processed.h5ad")
    slice2_path: str = str(_DATA_ROOT / "adata_118_processed.h5ad")
    preprocess: bool = True
    recipe: str = "generic"
    spatial_key: str = "spatial"
    counts_layer: str = "counts"
    label_key: str = "celltype"
    key_added: str = "align_spatial"
    min_genes: int = 10
    min_cells: int = 3
    feature_method: str = "hvg"
    n_top_genes: int = 2000
    run_pca: Optional[bool] = None
    pca_key: str = "X_pca"
    use_pca: bool = True
    dissimilarity: str = "cos"
    beta: float = 1.0
    lambdaVF: float = 1.0
    max_iter: int = 300
    K: int = 100
    use_downsampling: bool = False
    n_sampling: int = 10000
    sampling_method: str = "trn"
    device: Optional[str] = None
    cuda_visible_devices: Optional[str] = None
    outdir: Optional[str] = None
    alignment_kwargs: dict = field(default_factory=dict)


@dataclass
class TwoSliceResult:
    slice1_path: str = ""
    slice2_path: str = ""
    aligned_key: str = ""
    device: str = ""
    n_obs_slice1: int = 0
    n_obs_slice2: int = 0
    outdir: str = ""


# ---------------------------------------------------------------------------
# Atomic helpers
# ---------------------------------------------------------------------------


def load_two_slices(
    slice1_path: PathLike,
    slice2_path: PathLike,
) -> tuple[Any, Any]:
    p1, p2 = Path(slice1_path), Path(slice2_path)
    if not p1.exists():
        raise FileNotFoundError(f"slice1 not found: {p1}")
    if not p2.exists():
        raise FileNotFoundError(f"slice2 not found: {p2}")
    return load_h5ad(p1), load_h5ad(p2)


def preprocess_slices(
    slices: Sequence[Any],
    *,
    recipe: str = "generic",
    spatial_key: str = "spatial",
    counts_layer: str = "counts",
    min_genes: int = 10,
    min_cells: int = 3,
    feature_method: str = "hvg",
    n_top_genes: int = 2000,
    run_pca: Optional[bool] = None,
    **kwargs: Any,
) -> Sequence[Any]:
    from spateo.preprocessing.protocol_pipeline import preprocess_spatial

    for adata in slices:
        call: dict[str, Any] = dict(
            recipe=recipe,
            spatial_key=spatial_key,
            counts_layer=counts_layer,
            min_genes=min_genes,
            min_cells=min_cells,
            feature_method=feature_method,
            n_top_genes=n_top_genes,
        )
        if run_pca is not None:
            call["run_pca"] = run_pca
        call.update(kwargs)
        preprocess_spatial(adata, **call)
    return slices


def add_group_pca(
    slices: Sequence[Any],
    *,
    pca_key: str = "X_pca",
) -> Sequence[Any]:
    import spateo as st

    st.align.group_pca(list(slices), pca_key=pca_key)
    return slices


def align_two_slices(
    slice1: Any,
    slice2: Any,
    *,
    spatial_key: str = "spatial",
    key_added: str = "align_spatial",
    device: Optional[str] = None,
    use_pca: bool = True,
    pca_key: str = "X_pca",
    dissimilarity: str = "cos",
    beta: float = 1,
    lambdaVF: float = 1,
    max_iter: int = 300,
    K: int = 100,
    verbose: bool = True,
    **kwargs: Any,
) -> tuple[Sequence[Any], Any]:
    import spateo as st

    device = device or infer_device()
    models = [slice1, slice2]
    if use_pca:
        add_group_pca(models, pca_key=pca_key)
        kwargs.setdefault("rep_layer", pca_key)
        kwargs.setdefault("rep_field", "obsm")
        kwargs.setdefault("dissimilarity", dissimilarity)

    return st.align.morpho_align(
        models=models,
        spatial_key=spatial_key,
        key_added=key_added,
        device=device,
        verbose=verbose,
        beta=beta,
        lambdaVF=lambdaVF,
        max_iter=max_iter,
        K=K,
        **kwargs,
    )


def align_two_slices_downsampled(
    slice1: Any,
    slice2: Any,
    *,
    n_sampling: int = 10000,
    sampling_method: str = "trn",
    spatial_key: str = "spatial",
    key_added: str = "align_spatial",
    device: Optional[str] = None,
    use_pca: bool = True,
    pca_key: str = "X_pca",
    dissimilarity: str = "cos",
    beta: float = 1,
    lambdaVF: float = 1,
    max_iter: int = 300,
    K: int = 100,
    verbose: bool = True,
    **kwargs: Any,
) -> tuple[Sequence[Any], Sequence[Any], Any, Any]:
    import spateo as st

    device = device or infer_device()
    models = [slice1, slice2]
    if use_pca:
        add_group_pca(models, pca_key=pca_key)
        kwargs.setdefault("rep_layer", pca_key)
        kwargs.setdefault("rep_field", "obsm")
        kwargs.setdefault("dissimilarity", dissimilarity)

    return st.align.morpho_align_ref(
        models=models,
        spatial_key=spatial_key,
        n_sampling=n_sampling,
        sampling_method=sampling_method,
        key_added=key_added,
        device=device,
        verbose=verbose,
        beta=beta,
        lambdaVF=lambdaVF,
        max_iter=max_iter,
        K=K,
        **kwargs,
    )


def plot_input_slices(
    slices: Sequence[Any],
    *,
    label_key: str = "celltype",
    spatial_key: str = "spatial",
    height: int = 4,
    **kwargs: Any,
) -> Any:
    import spateo as st

    return st.pl.slices_2d(
        slices=list(slices),
        label_key=label_key,
        spatial_key=spatial_key,
        height=height,
        center_coordinate=True,
        show_legend=True,
        **kwargs,
    )


def plot_aligned_overlay(
    aligned_slices: Sequence[Any],
    *,
    spatial_key: str = "align_spatial_nonrigid",
    height: int = 3,
    overlay_type: str = "backward",
    **kwargs: Any,
) -> Any:
    import spateo as st

    return st.pl.overlay_slices_2d(
        slices=list(aligned_slices),
        spatial_key=spatial_key,
        height=height,
        overlay_type=overlay_type,
        **kwargs,
    )


def save_optimization_animation(
    aligned_slices: Sequence[Any],
    filename: PathLike,
    *,
    spatial_key: str = "spatial",
    key_added: str = "align_spatial",
    iter_key_added: str = "iter_spatial",
    label_key: str = "celltype",
    fps: int = 10,
    stepsize: int = 10,
    **kwargs: Any,
) -> Any:
    import spateo as st

    return st.pl.optimization_animation(
        aligned_slices=list(aligned_slices),
        spatial_key=spatial_key,
        key_added=key_added,
        iter_key_added=iter_key_added,
        filename=str(filename),
        fps=fps,
        stepsize=stepsize,
        label_key=label_key,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def run_two_slice_pipeline(config: TwoSliceConfig) -> TwoSliceResult:
    device = infer_device(config.cuda_visible_devices) if config.device is None else config.device
    LOGGER.info("Two-slice alignment: device=%s", device)

    slice1, slice2 = load_two_slices(config.slice1_path, config.slice2_path)
    validate_spatial(slice1, config.spatial_key)
    validate_spatial(slice2, config.spatial_key)
    LOGGER.info("Loaded slices: n_obs=%d, %d", slice1.n_obs, slice2.n_obs)

    slices = [slice1, slice2]
    if config.preprocess:
        preprocess_slices(
            slices,
            recipe=config.recipe,
            spatial_key=config.spatial_key,
            counts_layer=config.counts_layer,
            min_genes=config.min_genes,
            min_cells=config.min_cells,
            feature_method=config.feature_method,
            n_top_genes=config.n_top_genes,
            run_pca=config.run_pca,
        )
    add_group_pca(slices, pca_key=config.pca_key)

    params = dict(config.alignment_kwargs)
    if config.use_downsampling:
        aligned_slices, aligned_refs, extra1, extra2 = align_two_slices_downsampled(
            slice1,
            slice2,
            n_sampling=config.n_sampling,
            sampling_method=config.sampling_method,
            spatial_key=config.spatial_key,
            key_added=config.key_added,
            device=device,
            use_pca=config.use_pca,
            pca_key=config.pca_key,
            dissimilarity=config.dissimilarity,
            beta=config.beta,
            lambdaVF=config.lambdaVF,
            max_iter=config.max_iter,
            K=config.K,
            **params,
        )
    else:
        aligned_slices, pis = align_two_slices(
            slice1,
            slice2,
            spatial_key=config.spatial_key,
            key_added=config.key_added,
            device=device,
            use_pca=config.use_pca,
            pca_key=config.pca_key,
            dissimilarity=config.dissimilarity,
            beta=config.beta,
            lambdaVF=config.lambdaVF,
            max_iter=config.max_iter,
            K=config.K,
            **params,
        )

    outdir = ""
    if config.outdir:
        out = Path(config.outdir)
        out.mkdir(parents=True, exist_ok=True)
        outdir = str(out)

    return TwoSliceResult(
        slice1_path=str(config.slice1_path),
        slice2_path=str(config.slice2_path),
        aligned_key=config.key_added,
        device=device,
        n_obs_slice1=slice1.n_obs,
        n_obs_slice2=slice2.n_obs,
        outdir=outdir,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Two-slice spatial alignment")
    p.add_argument("--slice1", type=str, default=None)
    p.add_argument("--slice2", type=str, default=None)
    p.add_argument("--no-preprocess", action="store_true")
    p.add_argument("--spatial-key", type=str, default="spatial")
    p.add_argument("--key-added", type=str, default="align_spatial")
    p.add_argument("--use-downsampling", action="store_true")
    p.add_argument("--n-sampling", type=int, default=10000)
    p.add_argument("--device", type=str, default=None)
    p.add_argument("--outdir", type=str, default=None)
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_arg_parser().parse_args(argv)
    configure_logging(args.verbose)

    cfg = TwoSliceConfig(
        preprocess=not args.no_preprocess,
        spatial_key=args.spatial_key,
        key_added=args.key_added,
        use_downsampling=args.use_downsampling,
        n_sampling=args.n_sampling,
        device=args.device,
        outdir=args.outdir,
    )
    if args.slice1:
        cfg.slice1_path = args.slice1
    if args.slice2:
        cfg.slice2_path = args.slice2

    result = run_two_slice_pipeline(cfg)
    print(json.dumps(asdict(result), indent=2, ensure_ascii=False))


__all__ = [
    "TwoSliceConfig",
    "TwoSliceResult",
    "load_two_slices",
    "preprocess_slices",
    "add_group_pca",
    "align_two_slices",
    "align_two_slices_downsampled",
    "plot_input_slices",
    "plot_aligned_overlay",
    "save_optimization_animation",
    "run_two_slice_pipeline",
    "build_arg_parser",
    "main",
]
