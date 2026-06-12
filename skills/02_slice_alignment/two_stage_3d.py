"""Two-stage 3D model alignment skill.

Source: Alignment of 3D models(skill_code).ipynb
Data:   ../../data/skills_data/3.3D model alignment/
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Union

import numpy as np

from ..shared.data_helpers import (
    LOGGER,
    PathLike,
    configure_logging,
    infer_device,
    load_h5ad,
    validate_spatial,
)

_DATA_ROOT = Path(__file__).resolve().parent.parent.parent / "data" / "skills_data" / "3.3D model alignment"

DEFAULT_TISSUE_COLORMAP: dict[str, str] = {
    "muscle": "#5A2686",
    "epidermis": "#8A2BE2",
    "CNS": "#FF69B4",
    "midgut": "#DC143C",
    "yolk": "#B22222",
    "PNS": "#FF1493",
    "mesoderm": "#C71585",
    "trachea": "#FFA500",
    "amnioserosa": "#FFD700",
    "unknown": "#ADD8E6",
    "salivary gland": "#7F90F0",
    "hemolymph": "#6A5ACD",
}


# ---------------------------------------------------------------------------
# Config / Result
# ---------------------------------------------------------------------------


@dataclass
class TwoStage3DConfig:
    stage1_path: str = str(_DATA_ROOT / "07.62h_E8-10h_b_new_final_modified.h5ad")
    stage2_path: str = str(_DATA_ROOT / "08.84h_E4-8h_g_new_final_modified.h5ad")
    set_type: Optional[str] = "UMI"
    source_spatial_key: str = "align_spatial"
    before_align_key: str = "before_3d_align_spatial"
    aligned_key: str = "3d_align_spatial"
    counts_layer: str = "counts_X"
    log_layer: str = "log1p_X"
    groupby: str = "Annotation_2_tissue"
    key_added: str = "tissue"
    obs_index_key: str = "obs_index"
    stage1_mc_scale_factor: float = 1.2
    stage2_mc_scale_factor: float = 1.4
    stage1_surface_scale: float = 1.02
    stage2_surface_scale: float = 1.08
    alpha: float = 0.6
    smooth: int = 8000
    rotate_angle: tuple[float, float, float] = (165, 0, -10)
    n_sampling: int = 2000
    sampling_method: str = "trn"
    rep_layer: str = "log1p_X"
    rep_field: str = "layer"
    device: str = "0"
    outdir: Optional[str] = None
    align_kwargs: dict = field(default_factory=dict)


@dataclass
class TwoStage3DResult:
    stage1_path: str = ""
    stage2_path: str = ""
    aligned_key: str = ""
    n_obs_stage1: int = 0
    n_obs_stage2: int = 0
    saved_models: list[str] = field(default_factory=list)
    saved_h5ad: list[str] = field(default_factory=list)
    outdir: str = ""


# ---------------------------------------------------------------------------
# Atomic helpers
# ---------------------------------------------------------------------------


def read_stage_adata(path: PathLike, *, set_type: Optional[str] = "UMI") -> Any:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"stage adata not found: {p}")
    adata = load_h5ad(p)
    if set_type is not None:
        adata.uns["__type"] = set_type
    return adata


def ensure_log1p_layer(
    adata: Any,
    *,
    counts_layer: str = "counts_X",
    log_layer: str = "log1p_X",
) -> Any:
    if log_layer not in adata.layers:
        if counts_layer not in adata.layers:
            raise KeyError(f"counts layer {counts_layer!r} not in adata.layers")
        adata.layers[log_layer] = np.log1p(adata.layers[counts_layer])
    return adata


def reconstruct_3d_model(
    adata: Any,
    *,
    spatial_key: str = "align_spatial",
    groupby: str = "Annotation_2_tissue",
    key_added: str = "tissue",
    colormap: Any = "rainbow",
    alpha: float = 0.6,
    mc_scale_factor: float = 1.2,
    smooth: int = 8000,
    scale_factor: float = 1.02,
    construct_mesh: bool = True,
    **kwargs: Any,
) -> dict[str, Any]:
    import spateo as st

    pc, plot_cmap = st.tdr.construct_pc(
        adata=adata.copy(),
        spatial_key=spatial_key,
        groupby=groupby,
        key_added=key_added,
        colormap=colormap,
    )
    mesh = None
    if construct_mesh:
        mesh, _, _ = st.tdr.construct_surface(
            pc=pc,
            key_added=key_added,
            alpha=alpha,
            cs_method="marching_cube",
            cs_args={"mc_scale_factor": mc_scale_factor},
            smooth=smooth,
            scale_factor=scale_factor,
        )
    return {"pc": pc, "mesh": mesh, "colormap": plot_cmap}


def save_models(models: Mapping[str, Any], outdir: PathLike) -> dict[str, str]:
    import spateo as st

    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)
    saved: dict[str, str] = {}
    for name, model in models.items():
        if model is None:
            continue
        filename = outdir_p / (name if Path(name).suffix else f"{name}.vtk")
        st.tdr.save_model(model=model, filename=str(filename))
        saved[name] = str(filename)
    return saved


def read_models(model_paths: Mapping[str, PathLike]) -> dict[str, Any]:
    import spateo as st

    return {name: st.tdr.read_model(str(Path(path))) for name, path in model_paths.items()}


def rotate_stage_models(
    pc: Any,
    mesh: Optional[Any] = None,
    *,
    angle: tuple[float, float, float] = (165, 0, -10),
    rotate_center: Optional[Sequence[float]] = None,
) -> dict[str, Any]:
    import spateo as st

    center = rotate_center if rotate_center is not None else pc.center
    rotated_pc = st.tdr.rotate_model(model=pc, angle=angle, rotate_center=center, inplace=False)
    rotated_mesh = None
    if mesh is not None:
        rotated_mesh = st.tdr.rotate_model(model=mesh, angle=angle, rotate_center=center, inplace=False)
    return {"pc": rotated_pc, "mesh": rotated_mesh}


def prepare_stage_for_alignment(
    stage1_raw_adata: Any,
    stage1_rotated_pc: Any,
    stage2_adata: Any,
    *,
    obs_index_key: str = "obs_index",
    source_spatial_key: str = "align_spatial",
    before_align_key: str = "before_3d_align_spatial",
    counts_layer: str = "counts_X",
    log_layer: str = "log1p_X",
) -> tuple[Any, Any]:
    stage1_adata = stage1_raw_adata[np.asarray(stage1_rotated_pc.point_data[obs_index_key]), :].copy()
    ensure_log1p_layer(stage1_adata, counts_layer=counts_layer, log_layer=log_layer)
    stage1_adata.obsm[before_align_key] = np.asarray(stage1_rotated_pc.points)

    stage2_adata = stage2_adata.copy()
    ensure_log1p_layer(stage2_adata, counts_layer=counts_layer, log_layer=log_layer)
    stage2_adata.obsm[before_align_key] = stage2_adata.obsm[source_spatial_key].copy()
    return stage1_adata, stage2_adata


def align_two_stage_3d_models(
    stage1_adata: Any,
    stage2_adata: Any,
    *,
    n_sampling: int = 2000,
    sampling_method: str = "trn",
    rep_layer: str = "log1p_X",
    rep_field: str = "layer",
    spatial_key: str = "before_3d_align_spatial",
    key_added: str = "3d_align_spatial",
    device: Union[str, int] = "0",
    **kwargs: Any,
) -> dict[str, Any]:
    import spateo as st

    align_samples, align_samples_ref, extra1, extra2 = st.align.morpho_align_ref(
        models=[stage1_adata, stage2_adata],
        models_ref=None,
        n_sampling=n_sampling,
        sampling_method=sampling_method,
        rep_layer=rep_layer,
        rep_field=rep_field,
        spatial_key=spatial_key,
        key_added=key_added,
        device=device,
        **kwargs,
    )
    return {
        "stage1_aligned": align_samples[0].copy(),
        "stage2_aligned": align_samples[1].copy(),
        "stage1_aligned_ref": align_samples_ref[0].copy(),
        "stage2_aligned_ref": align_samples_ref[1].copy(),
        "raw": (align_samples, align_samples_ref, extra1, extra2),
    }


def build_aligned_point_clouds(
    alignment_result: Mapping[str, Any],
    *,
    before_align_key: str = "before_3d_align_spatial",
    aligned_key: str = "3d_align_spatial",
    groupby: str = "Annotation_2_tissue",
    key_added: str = "tissue",
    colormap: Any = DEFAULT_TISSUE_COLORMAP,
) -> dict[str, Any]:
    import spateo as st

    specs = {
        "stage1_raw_ref_pc": ("stage1_aligned_ref", before_align_key),
        "stage1_aligned_ref_pc": ("stage1_aligned_ref", aligned_key),
        "stage1_raw_pc": ("stage1_aligned", before_align_key),
        "stage1_aligned_pc": ("stage1_aligned", aligned_key),
        "stage2_raw_ref_pc": ("stage2_aligned_ref", before_align_key),
        "stage2_aligned_ref_pc": ("stage2_aligned_ref", aligned_key),
        "stage2_raw_pc": ("stage2_aligned", before_align_key),
        "stage2_aligned_pc": ("stage2_aligned", aligned_key),
    }
    outputs: dict[str, Any] = {}
    for out_name, (adata_key, sp_key) in specs.items():
        pc, _ = st.tdr.construct_pc(
            adata=alignment_result[adata_key].copy(),
            spatial_key=sp_key,
            groupby=groupby,
            key_added=key_added,
            colormap=colormap,
        )
        outputs[out_name] = pc
    return outputs


def save_aligned_adatas(
    alignment_result: Mapping[str, Any],
    outdir: PathLike,
    *,
    remove_iter_spatial_for_stage2: bool = False,
) -> dict[str, str]:
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)
    outputs = {
        "stage1_aligned": outdir_p / "stage1_3d_aligned.h5ad",
        "stage2_aligned": outdir_p / "stage2_3d_aligned.h5ad",
        "stage1_aligned_ref": outdir_p / "stage1_3d_aligned_ref.h5ad",
        "stage2_aligned_ref": outdir_p / "stage2_3d_aligned_ref.h5ad",
    }
    saved: dict[str, str] = {}
    for key, path in outputs.items():
        adata = alignment_result[key].copy()
        if remove_iter_spatial_for_stage2 and key == "stage2_aligned" and "iter_spatial" in adata.uns:
            del adata.uns["iter_spatial"]
        adata.write_h5ad(str(path))
        saved[key] = str(path)
    return saved


def plot_aligned_pair(
    pc1: Any,
    pc2: Any,
    *,
    key: str = "groups",
    colormap: Sequence[str] = ("#DC143C", "#0000FF"),
    **kwargs: Any,
) -> Any:
    import spateo as st

    pair = st.tdr.collect_models([pc1.copy(), pc2.copy()])
    defaults = dict(
        model_style="points", model_size=2, off_screen=True,
        show_legend=False, jupyter="static", show_axes=True,
    )
    defaults.update(kwargs)
    return st.pl.three_d_plot(model=pair, key=key, colormap=list(colormap), **defaults)


def label_and_plot_submodel_alignment(
    pcs: Sequence[Any],
    *,
    target_label: str = "amnioserosa",
    tissue_key: str = "tissue",
    output_key: str = "check_alignment",
    target_color: str = "#5A2686",
    other_color: str = "gainsboro",
    **plot_kwargs: Any,
) -> Any:
    import spateo as st

    aligned_pcs_v = []
    for pc in [p.copy() for p in pcs]:
        labels = np.asarray(pc.point_data[tissue_key]).copy()
        labels[labels != target_label] = "other"
        st.tdr.add_model_labels(
            model=pc,
            labels=labels,
            key_added=output_key,
            where="point_data",
            inplace=True,
            colormap={target_label: target_color, "other": other_color},
            alphamap={target_label: 1.0, "other": 0.5},
        )
        aligned_pcs_v.append(pc)
    defaults = dict(
        model_style="points", model_size=3, off_screen=False,
        shape=(1, len(aligned_pcs_v)), jupyter="static",
    )
    defaults.update(plot_kwargs)
    return st.pl.three_d_multi_plot(
        model=st.tdr.collect_models(aligned_pcs_v), key=output_key, **defaults,
    )


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def run_two_stage_3d_pipeline(config: TwoStage3DConfig) -> TwoStage3DResult:
    LOGGER.info("Two-stage 3D alignment: device=%s", config.device)

    stage1_raw = read_stage_adata(config.stage1_path, set_type=config.set_type)
    stage2_raw = read_stage_adata(config.stage2_path, set_type=config.set_type)
    validate_spatial(stage1_raw, config.source_spatial_key, require_3d=True)
    validate_spatial(stage2_raw, config.source_spatial_key, require_3d=True)
    LOGGER.info(
        "Loaded stages: n_obs=%d, %d", stage1_raw.n_obs, stage2_raw.n_obs,
    )

    stage1_models = reconstruct_3d_model(
        stage1_raw,
        spatial_key=config.source_spatial_key,
        groupby=config.groupby,
        key_added=config.key_added,
        mc_scale_factor=config.stage1_mc_scale_factor,
        scale_factor=config.stage1_surface_scale,
        alpha=config.alpha,
        smooth=config.smooth,
    )
    stage2_models = reconstruct_3d_model(
        stage2_raw,
        spatial_key=config.source_spatial_key,
        groupby=config.groupby,
        key_added=config.key_added,
        mc_scale_factor=config.stage2_mc_scale_factor,
        scale_factor=config.stage2_surface_scale,
        alpha=config.alpha,
        smooth=config.smooth,
    )

    rotated = rotate_stage_models(
        stage1_models["pc"], stage1_models["mesh"], angle=config.rotate_angle,
    )

    stage1_ready, stage2_ready = prepare_stage_for_alignment(
        stage1_raw,
        rotated["pc"],
        stage2_raw,
        obs_index_key=config.obs_index_key,
        source_spatial_key=config.source_spatial_key,
        before_align_key=config.before_align_key,
        counts_layer=config.counts_layer,
        log_layer=config.log_layer,
    )

    alignment = align_two_stage_3d_models(
        stage1_ready,
        stage2_ready,
        n_sampling=config.n_sampling,
        sampling_method=config.sampling_method,
        rep_layer=config.rep_layer,
        rep_field=config.rep_field,
        spatial_key=config.before_align_key,
        key_added=config.aligned_key,
        device=config.device,
        **dict(config.align_kwargs),
    )

    point_clouds = build_aligned_point_clouds(
        alignment,
        before_align_key=config.before_align_key,
        aligned_key=config.aligned_key,
        groupby=config.groupby,
        key_added=config.key_added,
    )

    saved_models: list[str] = []
    saved_h5ad: list[str] = []
    outdir = ""
    if config.outdir:
        out = Path(config.outdir)
        out.mkdir(parents=True, exist_ok=True)
        outdir = str(out)

        before = save_models(
            {
                "stage1_raw_pc_model.vtk": stage1_models["pc"],
                "stage1_raw_mesh_model.vtk": stage1_models["mesh"],
                "stage2_raw_pc_model.vtk": stage2_models["pc"],
                "stage2_raw_mesh_model.vtk": stage2_models["mesh"],
                "stage1_rotate_pc.vtk": rotated["pc"],
                "stage1_rotate_mesh.vtk": rotated["mesh"],
            },
            out / "models_before_align",
        )
        saved_models.extend(before.values())

        h5ad_paths = save_aligned_adatas(alignment, out / "aligned_h5ad")
        saved_h5ad.extend(h5ad_paths.values())

        aligned_models = save_models(point_clouds, out / "models_aligned")
        saved_models.extend(aligned_models.values())

    return TwoStage3DResult(
        stage1_path=str(config.stage1_path),
        stage2_path=str(config.stage2_path),
        aligned_key=config.aligned_key,
        n_obs_stage1=stage1_ready.n_obs,
        n_obs_stage2=stage2_ready.n_obs,
        saved_models=saved_models,
        saved_h5ad=saved_h5ad,
        outdir=outdir,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Two-stage 3D model alignment")
    p.add_argument("--stage1", type=str, default=None)
    p.add_argument("--stage2", type=str, default=None)
    p.add_argument("--source-spatial-key", type=str, default="align_spatial")
    p.add_argument("--aligned-key", type=str, default="3d_align_spatial")
    p.add_argument("--groupby", type=str, default="Annotation_2_tissue")
    p.add_argument("--n-sampling", type=int, default=2000)
    p.add_argument("--device", type=str, default="0")
    p.add_argument("--rotate-angle", type=float, nargs=3, default=[165, 0, -10])
    p.add_argument("--outdir", type=str, default=None)
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_arg_parser().parse_args(argv)
    configure_logging(args.verbose)

    cfg = TwoStage3DConfig(
        source_spatial_key=args.source_spatial_key,
        aligned_key=args.aligned_key,
        groupby=args.groupby,
        n_sampling=args.n_sampling,
        device=args.device,
        rotate_angle=tuple(args.rotate_angle),  # type: ignore[arg-type]
        outdir=args.outdir,
    )
    if args.stage1:
        cfg.stage1_path = args.stage1
    if args.stage2:
        cfg.stage2_path = args.stage2

    result = run_two_stage_3d_pipeline(cfg)
    print(json.dumps(asdict(result), indent=2, ensure_ascii=False))


__all__ = [
    "DEFAULT_TISSUE_COLORMAP",
    "TwoStage3DConfig",
    "TwoStage3DResult",
    "read_stage_adata",
    "ensure_log1p_layer",
    "reconstruct_3d_model",
    "save_models",
    "read_models",
    "rotate_stage_models",
    "prepare_stage_for_alignment",
    "align_two_stage_3d_models",
    "build_aligned_point_clouds",
    "save_aligned_adatas",
    "plot_aligned_pair",
    "label_and_plot_submodel_alignment",
    "run_two_stage_3d_pipeline",
    "build_arg_parser",
    "main",
]
