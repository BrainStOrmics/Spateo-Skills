"""3D reconstruction skill — point cloud, surface mesh, cell mesh, voxel, and subtype models."""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

warnings.filterwarnings("ignore")

LOGGER = logging.getLogger("spateo.reconstruction")

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# Lazy imports
# ---------------------------------------------------------------------------


def _import_spateo() -> Any:
    try:
        import spateo as st  # type: ignore
        return st
    except ImportError as exc:
        raise ImportError(
            "spateo is required for 3D reconstruction. "
            "Install from envs/requirements.txt."
        ) from exc


# ---------------------------------------------------------------------------
# Config / Result
# ---------------------------------------------------------------------------


@dataclass
class ReconstructionConfig:
    input_path: str
    out_dir: str
    prefix: str = "tdr_model"
    spatial_key: str = "3d_align_spatial"
    groupby: Optional[str] = None
    area_key: Optional[str] = None
    type_uns_key: str = "__type"
    type_uns_value: str = "UMI"

    models: List[str] = field(default_factory=lambda: ["point-cloud", "surface"])
    subtypes: List[str] = field(default_factory=list)

    point_label_key: str = "region"
    colormap: str = "rainbow"
    point_model_style: str = "points"

    surface_alpha: float = 0.6
    surface_method: str = "marching_cube"
    surface_mc_scale_factor: float = 0.8
    surface_smooth: int = 5000
    surface_scale_factor: float = 1.08

    cell_geometry: str = "sphere"
    cell_radius_factor: float = 0.3
    cell_radius_key: str = "cell_radius"

    voxel_label: str = "voxel"
    voxel_color: str = "gainsboro"
    voxel_smooth: int = 500

    subtype_alpha: float = 0.6
    subtype_color: str = "purple"
    subtype_mc_scale_factor: float = 0.9
    subtype_smooth: int = 6000
    subtype_scale_factor: float = 1.0

    save_models: bool = True
    save_plots: bool = False
    plot_format: str = "pdf"
    jupyter: Union[str, bool] = "static"
    window_size: Tuple[int, int] = (400, 400)
    show_axes: bool = True
    camera_position: Optional[Any] = None
    start_xvfb: bool = False
    off_screen: bool = False
    overwrite: bool = False


@dataclass
class ReconstructionResult:
    output_dir: str
    point_cloud_model: Optional[str] = None
    surface_model: Optional[str] = None
    cell_mesh_model: Optional[str] = None
    voxel_model: Optional[str] = None
    subtype_models: Dict[str, Dict[str, str]] = field(default_factory=dict)
    figures: Dict[str, str] = field(default_factory=dict)
    metadata_json: Optional[str] = None
    notes: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def safe_name(value: Any, max_len: int = 120) -> str:
    text = re.sub(r"\s+", "_", str(value).strip())
    text = re.sub(r"[^A-Za-z0-9_.\-]+", "_", text).strip("._-") or "unnamed"
    return text[:max_len]


def ensure_dir(path: PathLike) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def check_output_path(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists and overwrite=False: {path}")


def parse_jupyter_value(value: str) -> Union[str, bool]:
    low = value.lower()
    if low in {"true", "1", "yes"}:
        return True
    if low in {"false", "0", "no", "none"}:
        return False
    return value


def parse_camera_position(camera_json: Optional[str]) -> Optional[Any]:
    if not camera_json:
        return None
    try:
        return json.loads(camera_json)
    except json.JSONDecodeError as exc:
        raise ValueError(
            "--camera-position-json must be valid JSON, e.g. '[[x,y,z],[x,y,z],[x,y,z]]'."
        ) from exc


def set_cuda_visible_devices(cuda_visible_devices: Optional[str]) -> None:
    if cuda_visible_devices is not None:
        os.environ["CUDA_VISIBLE_DEVICES"] = cuda_visible_devices
        LOGGER.info("Set CUDA_VISIBLE_DEVICES=%s", cuda_visible_devices)


def maybe_start_xvfb(start_xvfb: bool, off_screen: bool) -> None:
    if not (start_xvfb or off_screen):
        return
    try:
        import pyvista as pv  # type: ignore
        if off_screen:
            pv.OFF_SCREEN = True
        if start_xvfb:
            pv.start_xvfb()
    except Exception as exc:
        LOGGER.warning("Could not configure PyVista/Xvfb: %s", exc)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_adata_for_reconstruction(
    adata: Any,
    spatial_key: str,
    groupby: Optional[str],
    area_key: Optional[str] = None,
) -> None:
    if spatial_key not in adata.obsm:
        raise KeyError(
            f"Spatial key {spatial_key!r} not found in adata.obsm. "
            f"Available: {list(adata.obsm.keys())}"
        )
    coords = np.asarray(adata.obsm[spatial_key])
    if coords.ndim != 2 or coords.shape[1] < 3:
        raise ValueError(
            f"adata.obsm[{spatial_key!r}] must have shape (n, >=3) for 3D. Got {coords.shape}."
        )
    if not np.isfinite(coords[:, :3]).all():
        raise ValueError(f"adata.obsm[{spatial_key!r}] contains NaN or Inf values.")
    if groupby is not None and groupby not in adata.obs:
        raise KeyError(f"Group column {groupby!r} not in adata.obs.")
    if area_key is not None and area_key not in adata.obs:
        raise KeyError(f"Area column {area_key!r} not in adata.obs.")


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_adata(input_path: PathLike) -> Any:
    path = Path(input_path)
    if not path.exists():
        raise FileNotFoundError(f"Input h5ad not found: {path}")
    st = _import_spateo()
    LOGGER.info("Reading AnnData: %s", path)
    adata = st.read_h5ad(str(path))
    adata.uns["__type"] = "UMI"
    return adata


def ensure_group_column(adata: Any, groupby: Optional[str]) -> str:
    if groupby is not None:
        return groupby
    fallback = "__tdr_all_points"
    adata.obs[fallback] = "all"
    return fallback


# ---------------------------------------------------------------------------
# Model construction — atomic functions
# ---------------------------------------------------------------------------


def construct_point_cloud(
    adata: Any,
    spatial_key: str,
    groupby: Optional[str],
    key_added: str = "region",
    colormap: str = "rainbow",
    **kwargs: Any,
) -> Tuple[Any, Any, str]:
    st = _import_spateo()
    group_key = ensure_group_column(adata, groupby)
    LOGGER.info("Constructing point cloud: obsm[%s] grouped by obs[%s].", spatial_key, group_key)
    pc, plot_cmap = st.tdr.construct_pc(
        adata=adata.copy(),
        spatial_key=spatial_key,
        groupby=group_key,
        key_added=key_added,
        colormap=colormap,
        **kwargs,
    )
    return pc, plot_cmap, group_key


def construct_surface_mesh(
    point_cloud: Any,
    key_added: str = "region",
    alpha: float = 0.6,
    cs_method: str = "marching_cube",
    mc_scale_factor: float = 0.8,
    smooth: int = 5000,
    scale_factor: float = 1.08,
    label: Optional[str] = None,
    color: Optional[str] = None,
    **kwargs: Any,
) -> Tuple[Any, Any, Any]:
    st = _import_spateo()
    LOGGER.info("Constructing surface mesh.")
    call_kwargs: Dict[str, Any] = dict(
        pc=point_cloud,
        key_added=key_added,
        alpha=alpha,
        cs_method=cs_method,
        cs_args={"mc_scale_factor": mc_scale_factor},
        smooth=smooth,
        scale_factor=scale_factor,
        **kwargs,
    )
    if label is not None:
        call_kwargs["label"] = label
    if color is not None:
        call_kwargs["color"] = color
    return st.tdr.construct_surface(**call_kwargs)


def get_area_values_for_point_cloud(adata: Any, point_cloud: Any, area_key: str) -> np.ndarray:
    if "obs_index" not in point_cloud.point_data:
        raise KeyError("Point cloud missing point_data['obs_index']; cannot map area values.")
    obs_index = [
        v.decode("utf-8") if isinstance(v, bytes) else str(v)
        for v in point_cloud.point_data["obs_index"].tolist()
    ]
    obs_names = set(map(str, adata.obs_names))
    if all(idx in obs_names for idx in obs_index):
        area = adata.obs.loc[obs_index, area_key].astype(float).to_numpy()
    else:
        try:
            area = adata.obs.iloc[[int(i) for i in obs_index]][area_key].astype(float).to_numpy()
        except ValueError as exc:
            raise KeyError("Could not map obs_index to adata.obs.") from exc
    if not np.isfinite(area).all():
        raise ValueError(f"adata.obs[{area_key!r}] contains NaN or Inf.")
    if (area < 0).any():
        raise ValueError(f"adata.obs[{area_key!r}] contains negative values.")
    return area


def construct_cell_mesh(
    adata: Any,
    point_cloud: Any,
    area_key: str,
    key_added: str = "cell_radius",
    geometry: str = "sphere",
    factor: float = 0.3,
    colormap: str = "hot_r",
    **kwargs: Any,
) -> Any:
    st = _import_spateo()
    LOGGER.info("Constructing cell mesh from area column %r.", area_key)
    area = get_area_values_for_point_cloud(adata, point_cloud, area_key)
    cell_radius = np.sqrt(area)
    st.tdr.add_model_labels(
        model=point_cloud, labels=cell_radius,
        key_added=key_added, where="point_data", colormap=colormap, inplace=True,
    )
    return st.tdr.construct_cells(
        pc=point_cloud, cell_size=point_cloud.point_data[key_added],
        geometry=geometry, factor=factor, **kwargs,
    )


def construct_voxel_model(
    mesh: Any,
    key_added: str = "region",
    label: str = "voxel",
    color: str = "gainsboro",
    smooth: int = 500,
    **kwargs: Any,
) -> Tuple[Any, Any]:
    st = _import_spateo()
    LOGGER.info("Voxelizing surface mesh.")
    return st.tdr.voxelize_mesh(
        mesh=mesh, voxel_pc=None, key_added=key_added,
        label=label, color=color, smooth=smooth, **kwargs,
    )


def pick_group_point_cloud(point_cloud: Any, key: str, group: str, **kwargs: Any) -> Any:
    st = _import_spateo()
    LOGGER.info("Selecting subtype/group: %s", group)
    return st.tdr.three_d_pick(model=point_cloud, key=key, picked_groups=group, **kwargs)[0]


def collect_models(models: Sequence[Any]) -> Any:
    st = _import_spateo()
    return st.tdr.collect_models(list(models))


# ---------------------------------------------------------------------------
# Saving and plotting
# ---------------------------------------------------------------------------


def save_model(model: Any, filename: PathLike, overwrite: bool = False) -> str:
    st = _import_spateo()
    filename = Path(filename)
    ensure_dir(filename.parent)
    check_output_path(filename, overwrite)
    LOGGER.info("Saving model: %s", filename)
    st.tdr.save_model(model=model, filename=str(filename))
    return str(filename)


def plot_model(
    model: Any,
    filename: PathLike,
    key: str = "region",
    model_style: Union[str, Sequence[str]] = "points",
    model_size: Optional[int] = None,
    colormap: Optional[Any] = None,
    show_axes: bool = True,
    jupyter: Union[str, bool] = "static",
    window_size: Tuple[int, int] = (400, 400),
    cpo: Optional[Any] = None,
    overwrite: bool = False,
    **kwargs: Any,
) -> str:
    st = _import_spateo()
    filename = Path(filename)
    ensure_dir(filename.parent)
    check_output_path(filename, overwrite)
    LOGGER.info("Rendering figure: %s", filename)
    call_kwargs: Dict[str, Any] = dict(
        model=model, key=key, model_style=model_style,
        show_axes=show_axes, jupyter=jupyter, window_size=window_size,
        cpo=cpo, filename=str(filename), **kwargs,
    )
    if model_size is not None:
        call_kwargs["model_size"] = model_size
    if colormap is not None:
        call_kwargs["colormap"] = colormap
    st.pl.three_d_plot(**call_kwargs)
    return str(filename)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def run_reconstruction_pipeline(config: ReconstructionConfig) -> ReconstructionResult:
    ensure_dir(config.out_dir)
    models_dir = ensure_dir(Path(config.out_dir) / "models")
    figures_dir = ensure_dir(Path(config.out_dir) / "figures")
    result = ReconstructionResult(output_dir=str(config.out_dir))

    maybe_start_xvfb(config.start_xvfb, config.off_screen)

    adata = load_adata(config.input_path)
    if config.type_uns_key:
        adata.uns[config.type_uns_key] = config.type_uns_value

    validate_adata_for_reconstruction(
        adata, config.spatial_key, config.groupby,
        area_key=config.area_key if "cell-mesh" in config.models else None,
    )

    pc, plot_cmap, actual_group_key = construct_point_cloud(
        adata, config.spatial_key, config.groupby,
        key_added=config.point_label_key, colormap=config.colormap,
    )

    prefix = safe_name(config.prefix)

    if config.save_models or "point-cloud" in config.models:
        result.point_cloud_model = save_model(
            pc, models_dir / f"{prefix}_point_cloud.vtk", overwrite=config.overwrite,
        )
    if config.save_plots:
        result.figures["point_cloud"] = plot_model(
            model=pc, key=config.point_label_key,
            model_style=config.point_model_style, colormap=plot_cmap,
            show_axes=config.show_axes, jupyter=config.jupyter,
            window_size=config.window_size, cpo=config.camera_position,
            filename=figures_dir / f"{prefix}_point_cloud.{config.plot_format}",
            overwrite=config.overwrite,
        )

    surface_mesh = None
    if any(m in config.models for m in ["surface", "voxel", "subtypes"]):
        surface_mesh, _, _ = construct_surface_mesh(
            pc, key_added=config.point_label_key, alpha=config.surface_alpha,
            cs_method=config.surface_method, mc_scale_factor=config.surface_mc_scale_factor,
            smooth=config.surface_smooth, scale_factor=config.surface_scale_factor,
        )
        if config.save_models or "surface" in config.models:
            result.surface_model = save_model(
                surface_mesh, models_dir / f"{prefix}_surface_mesh.vtk",
                overwrite=config.overwrite,
            )
        if config.save_plots:
            combined = collect_models([surface_mesh, pc])
            result.figures["surface_with_point_cloud"] = plot_model(
                model=combined, key=config.point_label_key,
                model_style=["surface", "points"], colormap=plot_cmap,
                show_axes=config.show_axes, jupyter=config.jupyter,
                window_size=config.window_size, cpo=config.camera_position,
                filename=figures_dir / f"{prefix}_surface_with_point_cloud.{config.plot_format}",
                overwrite=config.overwrite,
            )

    if "cell-mesh" in config.models:
        if not config.area_key:
            result.notes.append("Skipped cell-mesh: --area-key not provided.")
            LOGGER.warning("Skipping cell-mesh: area_key required.")
        else:
            cell_mesh = construct_cell_mesh(
                adata, pc, config.area_key,
                key_added=config.cell_radius_key, geometry=config.cell_geometry,
                factor=config.cell_radius_factor,
            )
            if config.save_models:
                result.cell_mesh_model = save_model(
                    cell_mesh, models_dir / f"{prefix}_cell_mesh.vtk",
                    overwrite=config.overwrite,
                )
            if config.save_plots:
                result.figures["cell_mesh"] = plot_model(
                    model=cell_mesh, key=config.point_label_key,
                    show_axes=config.show_axes, jupyter=config.jupyter,
                    window_size=config.window_size, cpo=config.camera_position,
                    filename=figures_dir / f"{prefix}_cell_mesh.{config.plot_format}",
                    overwrite=config.overwrite,
                )

    if "voxel" in config.models:
        if surface_mesh is None:
            raise RuntimeError("Surface mesh required for voxelization.")
        voxel_model, _ = construct_voxel_model(
            surface_mesh, key_added=config.point_label_key,
            label=config.voxel_label, color=config.voxel_color, smooth=config.voxel_smooth,
        )
        if config.save_models:
            result.voxel_model = save_model(
                voxel_model, models_dir / f"{prefix}_voxel.vtk", overwrite=config.overwrite,
            )
        if config.save_plots:
            result.figures["voxel"] = plot_model(
                model=voxel_model, key=config.point_label_key,
                show_axes=config.show_axes, jupyter=config.jupyter,
                window_size=config.window_size, cpo=config.camera_position,
                filename=figures_dir / f"{prefix}_voxel.{config.plot_format}",
                overwrite=config.overwrite,
            )

    if "subtypes" in config.models:
        if not config.subtypes:
            result.notes.append("Subtypes requested but none provided.")
            LOGGER.warning("Skipping subtypes: no subtypes listed.")
        elif surface_mesh is None:
            raise RuntimeError("Global surface mesh required for subtype reconstruction.")
        else:
            for subtype in config.subtypes:
                subtype_name = safe_name(subtype)
                sub_result: Dict[str, str] = {}
                subtype_pc = pick_group_point_cloud(pc, config.point_label_key, subtype)
                subtype_mesh, subtype_surface_pc, _ = construct_surface_mesh(
                    subtype_pc, key_added=config.point_label_key,
                    label=subtype, color=config.subtype_color,
                    alpha=config.subtype_alpha, cs_method=config.surface_method,
                    mc_scale_factor=config.subtype_mc_scale_factor,
                    smooth=config.subtype_smooth, scale_factor=config.subtype_scale_factor,
                )
                if config.save_models:
                    sub_result["point_cloud"] = save_model(
                        subtype_surface_pc,
                        models_dir / f"{prefix}_{subtype_name}_point_cloud.vtk",
                        overwrite=config.overwrite,
                    )
                    sub_result["mesh"] = save_model(
                        subtype_mesh,
                        models_dir / f"{prefix}_{subtype_name}_mesh.vtk",
                        overwrite=config.overwrite,
                    )
                if config.save_plots:
                    combined = collect_models([surface_mesh, subtype_mesh, subtype_surface_pc])
                    sub_result["figure"] = plot_model(
                        model=combined, key=config.point_label_key,
                        model_style=["surface", "surface", "points"], model_size=3,
                        show_axes=config.show_axes, jupyter=config.jupyter,
                        window_size=config.window_size, cpo=config.camera_position,
                        filename=figures_dir / f"{prefix}_{subtype_name}_surface_mesh_point_cloud.{config.plot_format}",
                        overwrite=config.overwrite,
                    )
                result.subtype_models[subtype] = sub_result

    meta = {
        "config": _jsonable_config(config),
        "result": asdict(result),
        "adata_shape": list(adata.shape),
        "spatial_key": config.spatial_key,
        "groupby": actual_group_key,
        "models_requested": config.models,
    }
    meta_path = Path(config.out_dir) / f"{prefix}_reconstruction_metadata.json"
    check_output_path(meta_path, config.overwrite)
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    result.metadata_json = str(meta_path)

    LOGGER.info("3D reconstruction complete. Output: %s", config.out_dir)
    return result


def _jsonable_config(config: ReconstructionConfig) -> Dict[str, Any]:
    data = asdict(config)
    return {k: str(v) if isinstance(v, Path) else v for k, v in data.items()}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Build Spateo 3D reconstruction models from an AnnData with 3D coordinates.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", required=True, help="Input .h5ad file.")
    p.add_argument("--out-dir", required=True, help="Output directory.")
    p.add_argument("--prefix", default="tdr_model")
    p.add_argument("--spatial-key", default="3d_align_spatial")
    p.add_argument("--groupby", default=None)
    p.add_argument("--area-key", default=None)
    p.add_argument("--models", nargs="+", default=["point-cloud", "surface"],
                   choices=["point-cloud", "surface", "cell-mesh", "voxel", "subtypes"])
    p.add_argument("--subtypes", nargs="*", default=[])
    p.add_argument("--point-label-key", default="region")
    p.add_argument("--colormap", default="rainbow")
    p.add_argument("--surface-alpha", type=float, default=0.6)
    p.add_argument("--surface-method", default="marching_cube")
    p.add_argument("--surface-mc-scale-factor", type=float, default=0.8)
    p.add_argument("--surface-smooth", type=int, default=5000)
    p.add_argument("--surface-scale-factor", type=float, default=1.08)
    p.add_argument("--cell-geometry", default="sphere")
    p.add_argument("--cell-radius-factor", type=float, default=0.3)
    p.add_argument("--voxel-label", default="voxel")
    p.add_argument("--voxel-color", default="gainsboro")
    p.add_argument("--voxel-smooth", type=int, default=500)
    p.add_argument("--subtype-alpha", type=float, default=0.6)
    p.add_argument("--subtype-color", default="purple")
    p.add_argument("--subtype-mc-scale-factor", type=float, default=0.9)
    p.add_argument("--subtype-smooth", type=int, default=6000)
    p.add_argument("--subtype-scale-factor", type=float, default=1.0)
    p.add_argument("--save-models", action="store_true", default=True)
    p.add_argument("--no-save-models", action="store_false", dest="save_models")
    p.add_argument("--save-plots", action="store_true")
    p.add_argument("--plot-format", default="pdf")
    p.add_argument("--jupyter", default="static")
    p.add_argument("--window-size", nargs=2, type=int, default=[400, 400])
    p.add_argument("--camera-position-json", default=None)
    p.add_argument("--hide-axes", action="store_true")
    p.add_argument("--xvfb", action="store_true")
    p.add_argument("--off-screen", action="store_true")
    p.add_argument("--cuda-visible-devices", default=None)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    configure_logging(verbose=args.verbose)
    set_cuda_visible_devices(args.cuda_visible_devices)
    try:
        config = ReconstructionConfig(
            input_path=args.input,
            out_dir=args.out_dir,
            prefix=args.prefix,
            spatial_key=args.spatial_key,
            groupby=args.groupby,
            area_key=args.area_key,
            models=list(args.models),
            subtypes=list(args.subtypes),
            point_label_key=args.point_label_key,
            colormap=args.colormap,
            surface_alpha=args.surface_alpha,
            surface_method=args.surface_method,
            surface_mc_scale_factor=args.surface_mc_scale_factor,
            surface_smooth=args.surface_smooth,
            surface_scale_factor=args.surface_scale_factor,
            cell_geometry=args.cell_geometry,
            cell_radius_factor=args.cell_radius_factor,
            voxel_label=args.voxel_label,
            voxel_color=args.voxel_color,
            voxel_smooth=args.voxel_smooth,
            subtype_alpha=args.subtype_alpha,
            subtype_color=args.subtype_color,
            subtype_mc_scale_factor=args.subtype_mc_scale_factor,
            subtype_smooth=args.subtype_smooth,
            subtype_scale_factor=args.subtype_scale_factor,
            save_models=args.save_models,
            save_plots=args.save_plots,
            plot_format=args.plot_format,
            jupyter=parse_jupyter_value(args.jupyter),
            window_size=(args.window_size[0], args.window_size[1]),
            show_axes=not args.hide_axes,
            camera_position=parse_camera_position(args.camera_position_json),
            start_xvfb=args.xvfb,
            off_screen=args.off_screen,
            overwrite=args.overwrite,
        )
        result = run_reconstruction_pipeline(config)
    except Exception as exc:
        LOGGER.exception("3D reconstruction failed: %s", exc)
        return 1
    print(json.dumps(asdict(result), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "ReconstructionConfig",
    "ReconstructionResult",
    "build_arg_parser",
    "collect_models",
    "configure_logging",
    "construct_cell_mesh",
    "construct_point_cloud",
    "construct_surface_mesh",
    "construct_voxel_model",
    "load_adata",
    "main",
    "pick_group_point_cloud",
    "plot_model",
    "run_reconstruction_pipeline",
    "save_model",
    "validate_adata_for_reconstruction",
]
