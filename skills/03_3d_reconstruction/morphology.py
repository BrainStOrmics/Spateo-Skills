"""3D model morphology skill — KDE density and morphology metrics."""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

warnings.filterwarnings("ignore")

LOGGER = logging.getLogger("spateo.morphology")

PathLike = Union[str, Path]


DEFAULT_CPO: List[Tuple[float, float, float]] = [
    (1367.3561390400112, 2417.302906013576, 1707.0488857080545),
    (347.35097122128604, 662.2007202773002, 339.2838334086367),
    (-0.45946273061616877, -0.3642657132348304, 0.810064497022462),
]


# ---------------------------------------------------------------------------
# Lazy imports
# ---------------------------------------------------------------------------


def _import_spateo() -> Any:
    try:
        import spateo as st  # type: ignore
        return st
    except ImportError as exc:
        raise ImportError("spateo is required for morphology analysis.") from exc


# ---------------------------------------------------------------------------
# Config / Result
# ---------------------------------------------------------------------------


@dataclass
class MorphologyConfig:
    pc_model_path: str
    mesh_model_path: str
    output_dir: str = "./spateo_3d_model_morphology_output"
    bandwidth: float = 5.0
    kde_key: str = "cells_kde"
    kde_colormap: str = "hot_r"
    unit_scale: float = 0.001
    run_kde: bool = True
    run_morphology: bool = True
    save_kde_model: bool = True
    save_scaled_models: bool = True
    plot_kde: bool = True
    cpo: Optional[List[Tuple[float, float, float]]] = field(default_factory=lambda: list(DEFAULT_CPO))
    window_size: Tuple[int, int] = (400, 400)
    opacity: float = 0.6
    model_style: str = "points"
    jupyter: str = "static"
    overwrite: bool = False


@dataclass
class MorphologyResult:
    output_dir: str
    kde_model_path: Optional[str] = None
    kde_figure_path: Optional[str] = None
    scaled_pc_model_path: Optional[str] = None
    scaled_mesh_model_path: Optional[str] = None
    morphology_json_path: Optional[str] = None
    morphology_csv_path: Optional[str] = None
    morphology_txt_path: Optional[str] = None
    config_json_path: Optional[str] = None
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


def ensure_path(path: PathLike, label: str) -> Path:
    p = Path(path).expanduser().resolve()
    if not p.exists():
        raise FileNotFoundError(f"{label} not found: {p}")
    return p


def ensure_dir(path: PathLike) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def check_output_path(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists and overwrite=False: {path}")


def parse_window_size(value: Union[str, Sequence[int]]) -> Tuple[int, int]:
    if isinstance(value, str):
        parts = [x.strip() for x in value.replace("x", ",").split(",") if x.strip()]
        if len(parts) != 2:
            raise ValueError("window_size must be like '400,400' or '400x400'.")
        return int(parts[0]), int(parts[1])
    if len(value) != 2:
        raise ValueError("window_size must have exactly two integers.")
    return int(value[0]), int(value[1])


def load_cpo(cpo_json: Optional[str]) -> Optional[List[Tuple[float, float, float]]]:
    if cpo_json is None or cpo_json == "":
        return list(DEFAULT_CPO)
    if cpo_json.lower() in {"none", "null", "false"}:
        return None
    candidate = Path(cpo_json).expanduser()
    raw = candidate.read_text(encoding="utf-8") if candidate.exists() else cpo_json
    data = json.loads(raw)
    if len(data) != 3 or any(len(row) != 3 for row in data):
        raise ValueError("cpo must be a 3x3 JSON array.")
    return [(float(r[0]), float(r[1]), float(r[2])) for r in data]


def to_jsonable(obj: Any) -> Any:
    if obj is None or isinstance(obj, (str, int, float, bool)):
        return obj
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Mapping):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [to_jsonable(v) for v in obj]
    if hasattr(obj, "to_dict"):
        try:
            return to_jsonable(obj.to_dict())
        except Exception:
            pass
    return repr(obj)


# ---------------------------------------------------------------------------
# Core operations
# ---------------------------------------------------------------------------


def load_models(config: MorphologyConfig) -> Tuple[Any, Any]:
    st = _import_spateo()
    pc_path = ensure_path(config.pc_model_path, "Point-cloud model")
    mesh_path = ensure_path(config.mesh_model_path, "Mesh model")
    LOGGER.info("Loading point cloud: %s", pc_path)
    LOGGER.info("Loading mesh: %s", mesh_path)
    return st.tdr.read_model(str(pc_path)), st.tdr.read_model(str(mesh_path))


def calculate_point_cloud_kde(
    pc: Any,
    config: MorphologyConfig,
    output_dir: PathLike,
    **kwargs: Any,
) -> Tuple[Any, Optional[str], Optional[str]]:
    st = _import_spateo()
    outdir = ensure_dir(output_dir)
    figures_dir = ensure_dir(outdir / "figures")

    LOGGER.info("Calculating point-cloud KDE (bandwidth=%.2f).", config.bandwidth)
    st.tdr.pc_KDE(
        pc=pc, bandwidth=config.bandwidth, key_added=config.kde_key,
        colormap=config.kde_colormap, inplace=True, **kwargs,
    )

    kde_model_path: Optional[str] = None
    if config.save_kde_model:
        kde_model_path = str(outdir / f"pc_model_with_{config.kde_key}.vtk")
        st.tdr.save_model(pc, kde_model_path)

    kde_figure_path: Optional[str] = None
    if config.plot_kde:
        kde_figure_path = str(figures_dir / f"{config.kde_key}.pdf")
        st.pl.three_d_plot(
            model=pc, key=config.kde_key, colormap=config.kde_colormap,
            opacity=config.opacity, model_style=config.model_style,
            jupyter=config.jupyter, window_size=config.window_size,
            cpo=config.cpo, filename=kde_figure_path,
        )

    return pc, kde_model_path, kde_figure_path


def scale_model_coordinates(model: Any, unit_scale: float) -> Any:
    scaled = model.copy()
    scaled.points = np.asarray(scaled.points) * float(unit_scale)
    return scaled


def calculate_model_morphology(
    pc: Any,
    mesh: Any,
    config: MorphologyConfig,
    output_dir: PathLike,
    **kwargs: Any,
) -> Tuple[Any, Optional[str], Optional[str], str, Optional[str], str]:
    st = _import_spateo()
    outdir = ensure_dir(output_dir)

    LOGGER.info("Scaling coordinates (factor=%.4f) and computing morphology.", config.unit_scale)
    scaled_pc = scale_model_coordinates(pc, config.unit_scale)
    scaled_mesh = scale_model_coordinates(mesh, config.unit_scale)

    scaled_pc_path: Optional[str] = None
    scaled_mesh_path: Optional[str] = None
    if config.save_scaled_models:
        scaled_pc_path = str(outdir / "pc_model_scaled.vtk")
        scaled_mesh_path = str(outdir / "mesh_model_scaled.vtk")
        st.tdr.save_model(scaled_pc, scaled_pc_path)
        st.tdr.save_model(scaled_mesh, scaled_mesh_path)

    morphology = st.tdr.model_morphology(model=scaled_mesh, pc=scaled_pc, **kwargs)

    json_path = outdir / "model_morphology.json"
    txt_path = outdir / "model_morphology.txt"
    txt_path.write_text(repr(morphology), encoding="utf-8")
    json_path.write_text(
        json.dumps(to_jsonable(morphology), ensure_ascii=False, indent=2), encoding="utf-8",
    )

    csv_path: Optional[str] = None
    if hasattr(morphology, "to_csv"):
        csv_out = outdir / "model_morphology.csv"
        try:
            morphology.to_csv(csv_out)
            csv_path = str(csv_out)
        except Exception:
            pass

    return morphology, scaled_pc_path, scaled_mesh_path, str(json_path), csv_path, str(txt_path)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def run_morphology_pipeline(config: MorphologyConfig) -> MorphologyResult:
    outdir = ensure_dir(config.output_dir)
    pc, mesh = load_models(config)
    result = MorphologyResult(output_dir=str(outdir))

    config_path = outdir / "config.json"
    config_path.write_text(
        json.dumps(to_jsonable(asdict(config)), ensure_ascii=False, indent=2), encoding="utf-8",
    )
    result.config_json_path = str(config_path)

    if config.run_kde:
        pc, kde_model_path, kde_figure_path = calculate_point_cloud_kde(pc, config, outdir)
        result.kde_model_path = kde_model_path
        result.kde_figure_path = kde_figure_path

    if config.run_morphology:
        (
            morphology,
            scaled_pc_path,
            scaled_mesh_path,
            json_path,
            csv_path,
            txt_path,
        ) = calculate_model_morphology(pc, mesh, config, outdir)
        result.scaled_pc_model_path = scaled_pc_path
        result.scaled_mesh_model_path = scaled_mesh_path
        result.morphology_json_path = json_path
        result.morphology_csv_path = csv_path
        result.morphology_txt_path = txt_path

    LOGGER.info("Morphology analysis complete. Output: %s", outdir)
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Calculate KDE and morphology features for Spateo 3D models.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--pc-model", required=True, help="Input point-cloud model (.vtk).")
    p.add_argument("--mesh-model", required=True, help="Input mesh/surface model (.vtk).")
    p.add_argument("--output-dir", default="./spateo_3d_model_morphology_output")
    p.add_argument("--bandwidth", type=float, default=5.0)
    p.add_argument("--kde-key", default="cells_kde")
    p.add_argument("--kde-colormap", default="hot_r")
    p.add_argument("--unit-scale", type=float, default=0.001,
                   help="Coordinate multiplier (default 0.001 = /1000).")
    p.add_argument("--no-kde", action="store_true")
    p.add_argument("--no-morphology", action="store_true")
    p.add_argument("--no-save-kde-model", action="store_true")
    p.add_argument("--no-save-scaled-models", action="store_true")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("--cpo-json", default=None,
                   help="Camera position JSON or file. 'none' to disable.")
    p.add_argument("--window-size", default="400,400")
    p.add_argument("--opacity", type=float, default=0.6)
    p.add_argument("--model-style", default="points")
    p.add_argument("--jupyter", default="static")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    configure_logging(verbose=args.verbose)
    try:
        config = MorphologyConfig(
            pc_model_path=args.pc_model,
            mesh_model_path=args.mesh_model,
            output_dir=args.output_dir,
            bandwidth=args.bandwidth,
            kde_key=args.kde_key,
            kde_colormap=args.kde_colormap,
            unit_scale=args.unit_scale,
            run_kde=not args.no_kde,
            run_morphology=not args.no_morphology,
            save_kde_model=not args.no_save_kde_model,
            save_scaled_models=not args.no_save_scaled_models,
            plot_kde=not args.no_plot,
            cpo=load_cpo(args.cpo_json),
            window_size=parse_window_size(args.window_size),
            opacity=args.opacity,
            model_style=args.model_style,
            jupyter=args.jupyter,
            overwrite=args.overwrite,
        )
        result = run_morphology_pipeline(config)
    except Exception as exc:
        LOGGER.exception("Morphology pipeline failed: %s", exc)
        return 1
    print(json.dumps(asdict(result), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "MorphologyConfig",
    "MorphologyResult",
    "build_arg_parser",
    "calculate_model_morphology",
    "calculate_point_cloud_kde",
    "configure_logging",
    "load_models",
    "main",
    "run_morphology_pipeline",
    "scale_model_coordinates",
]
