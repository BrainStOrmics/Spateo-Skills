"""3D transcriptomics interpolation skill — vtk, gp, kernel, and deep methods."""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

warnings.filterwarnings("ignore")

LOGGER = logging.getLogger("spateo.interpolation")

PathLike = Union[str, Path]


# ---------------------------------------------------------------------------
# Lazy imports
# ---------------------------------------------------------------------------


def _import_spateo() -> Any:
    try:
        import spateo as st  # type: ignore
        return st
    except ImportError as exc:
        raise ImportError("spateo is required for interpolation.") from exc


def _import_spateo_plus() -> Any:
    try:
        import spateo_plus as sp  # type: ignore
        return sp
    except ImportError as exc:
        raise ImportError("spateo_plus required for --normalize / --log1p.") from exc


def _import_sparse() -> Any:
    from scipy import sparse  # type: ignore
    return sparse


# ---------------------------------------------------------------------------
# Config / Result
# ---------------------------------------------------------------------------


@dataclass
class InterpolationConfig:
    input_path: str
    out_dir: str
    prefix: str = "tdr_interpolation"
    spatial_key: str = "3d_align_spatial"
    groupby: Optional[str] = None
    genes: List[str] = field(default_factory=list)
    top_n_genes: int = 0

    normalize: bool = False
    log1p: bool = False
    make_dense: bool = True

    sample_n: Optional[int] = None
    sample_method: str = "random"
    random_seed: int = 0

    methods: List[str] = field(default_factory=lambda: ["vtk", "gp", "kernel"])

    target_points_path: Optional[str] = None
    target_model_path: Optional[str] = None
    target_point_columns: Tuple[str, str, str] = ("x", "y", "z")
    build_voxel_target: bool = True

    surface_alpha: float = 0.6
    surface_method: str = "marching_cube"
    surface_mc_scale_factor: float = 1.0
    surface_smooth: int = 5000
    surface_scale_factor: float = 1.0
    voxel_smooth: int = 50

    vtk_n_points: int = 5
    gp_device: str = "cpu"

    save_interpolated_h5ad: bool = True
    save_models: bool = True
    save_target_models: bool = True
    save_plots: bool = False
    save_voxel_slices: bool = False
    slice_axis: str = "x"
    n_slices: int = 15
    slice_source_method: str = "gp"

    colormap: str = "hot_r"
    opacity: float = 0.5
    model_style: str = "points"
    jupyter: Union[str, bool] = "static"
    window_size: Tuple[int, int] = (500, 500)
    camera_position: Optional[Any] = None
    start_xvfb: bool = False
    off_screen: bool = False
    overwrite: bool = False


@dataclass
class MethodResult:
    method: str
    interpolated_h5ad: Optional[str] = None
    point_cloud_model: Optional[str] = None
    plot: Optional[str] = None
    voxel_slice_plot: Optional[str] = None
    error: Optional[str] = None


@dataclass
class InterpolationResult:
    output_dir: str
    genes: List[str]
    target_points_shape: Tuple[int, int]
    target_point_source: str
    raw_point_cloud_model: Optional[str] = None
    surface_model: Optional[str] = None
    voxel_model: Optional[str] = None
    methods: Dict[str, MethodResult] = field(default_factory=dict)
    metadata_path: Optional[str] = None
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


def configure_runtime(config: InterpolationConfig) -> None:
    if config.off_screen:
        os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
    if config.start_xvfb:
        try:
            import pyvista as pv  # type: ignore
            pv.start_xvfb()
        except Exception as exc:
            LOGGER.warning("Could not start Xvfb: %s", exc)


def ensure_out_dir(path: PathLike, overwrite: bool) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    if not overwrite:
        non_hidden = [f for f in p.iterdir() if not f.name.startswith(".")]
        if non_hidden:
            raise FileExistsError(f"Output dir not empty: {p}. Use --overwrite.")
    return p


def parse_camera_position(value: Optional[str]) -> Optional[Any]:
    if value is None:
        return None
    possible_path = Path(value)
    if possible_path.exists():
        value = possible_path.read_text(encoding="utf-8")
    try:
        return json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError("Camera position must be JSON or a path to a JSON file.") from exc


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_spatial(adata: Any, spatial_key: str) -> np.ndarray:
    if spatial_key not in adata.obsm:
        raise KeyError(
            f"adata.obsm[{spatial_key!r}] not found. Available: {list(adata.obsm.keys())}"
        )
    coords = np.asarray(adata.obsm[spatial_key])
    if coords.ndim != 2 or coords.shape[1] < 3:
        raise ValueError(
            f"adata.obsm[{spatial_key!r}] must be (n, >=3). Got shape {coords.shape}."
        )
    if not np.isfinite(coords[:, :3]).all():
        raise ValueError(f"adata.obsm[{spatial_key!r}] contains NaN or Inf.")
    return coords[:, :3]


def validate_genes(adata: Any, genes: Sequence[str]) -> None:
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise KeyError(f"Genes not in adata.var_names: {missing[:20]}")


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------


def load_adata(input_path: PathLike) -> Any:
    path = Path(input_path)
    if not path.exists():
        raise FileNotFoundError(f"Input h5ad not found: {path}")
    st = _import_spateo()
    LOGGER.info("Reading AnnData: %s", path)
    return st.read_h5ad(str(path))


def normalize_and_log1p(adata: Any, normalize: bool, log1p: bool) -> Any:
    if not normalize and not log1p:
        return adata
    sp = _import_spateo_plus()
    adata.uns.setdefault("pp", {})
    adata.uns.setdefault("__type", "UMI")
    if normalize:
        LOGGER.info("Normalizing AnnData.")
        sp.pp.normalize(adata)
    if log1p:
        LOGGER.info("Log-transforming AnnData.")
        sp.pp.log1p(adata)
    return adata


def downsample_adata(
    adata: Any,
    spatial_key: str,
    n: Optional[int],
    method: str = "random",
    random_seed: int = 0,
) -> Any:
    if n is None or n <= 0 or n >= adata.n_obs:
        return adata
    np.random.seed(random_seed)
    try:
        from spateo_plus.tools import sample  # type: ignore
        LOGGER.info("Downsampling with spateo_plus: n=%d, method=%s", n, method)
        sampling = sample(
            arr=np.asarray(adata.obs_names), n=n, method=method,
            X=np.asarray(adata.obsm[spatial_key])[:, :3],
        )
        return adata[sampling, :].copy()
    except Exception as exc:
        LOGGER.warning("spateo_plus sampling failed (%s); falling back to numpy.", exc)
        idx = np.random.choice(np.arange(adata.n_obs), size=n, replace=False)
        return adata[idx, :].copy()


def make_adata_dense(adata: Any) -> Any:
    sparse = _import_sparse()
    adata_dense = adata.copy()
    if sparse.issparse(adata_dense.X):
        LOGGER.info("Converting sparse .X to dense.")
        adata_dense.X = adata_dense.X.toarray()
    return adata_dense


def select_genes(adata: Any, genes: Sequence[str], top_n_genes: int = 0) -> List[str]:
    if genes:
        validate_genes(adata, genes)
        return list(genes)
    if top_n_genes <= 0:
        raise ValueError("Provide --genes or set --top-n-genes > 0.")
    x = adata.X
    gene_mean = np.asarray(x.mean(axis=0)).ravel() if hasattr(x, "toarray") else np.asarray(x).mean(axis=0).ravel()
    order = np.argsort(gene_mean)[::-1][:top_n_genes]
    selected = [str(adata.var_names[i]) for i in order]
    LOGGER.info("Top %d genes: %s", len(selected), selected)
    return selected


def expression_vector(adata: Any, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(f"Gene {gene!r} not in adata.var_names.")
    x = adata[:, gene].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).ravel()


# ---------------------------------------------------------------------------
# Target points
# ---------------------------------------------------------------------------


def read_target_points(path: PathLike, columns: Tuple[str, str, str]) -> np.ndarray:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Target points file not found: {path}")
    suffix = path.suffix.lower()
    if suffix == ".npy":
        points = np.load(path)
    elif suffix == ".npz":
        data = np.load(path)
        key = "points" if "points" in data.files else data.files[0]
        points = data[key]
    elif suffix in {".csv", ".tsv", ".txt"}:
        import pandas as pd  # type: ignore
        sep = "\t" if suffix in {".tsv", ".txt"} else ","
        df = pd.read_csv(path, sep=sep)
        if all(c in df.columns for c in columns):
            points = df.loc[:, list(columns)].to_numpy()
        else:
            numeric = df.select_dtypes(include=[np.number])
            if numeric.shape[1] < 3:
                raise ValueError(f"Cannot find columns {columns} or 3 numeric columns in {path}.")
            points = numeric.iloc[:, :3].to_numpy()
    else:
        raise ValueError(f"Unsupported target point format: {suffix}")
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] < 3:
        raise ValueError(f"Target points must be (n, >=3). Got {points.shape}.")
    return points[:, :3]


def build_voxel_target(
    adata: Any,
    config: InterpolationConfig,
) -> Tuple[np.ndarray, Any, Any, Any]:
    st = _import_spateo()
    LOGGER.info("Building voxel target from AnnData.")
    pc, _ = st.tdr.construct_pc(
        adata=adata.copy(), spatial_key=config.spatial_key, groupby=config.groupby,
    )
    mesh, _, _ = st.tdr.construct_surface(
        pc=pc, alpha=config.surface_alpha, cs_method=config.surface_method,
        cs_args={"mc_scale_factor": config.surface_mc_scale_factor},
        smooth=config.surface_smooth, scale_factor=config.surface_scale_factor,
    )
    voxel, _ = st.tdr.voxelize_mesh(mesh=mesh, voxel_pc=pc, smooth=config.voxel_smooth)
    return np.asarray(voxel.points), pc, mesh, voxel


def get_target_points(
    config: InterpolationConfig, adata: Any,
) -> Tuple[np.ndarray, str, Optional[Any], Optional[Any], Optional[Any]]:
    st = _import_spateo()
    if config.target_points_path is not None:
        pts = read_target_points(config.target_points_path, config.target_point_columns)
        return pts, "target_points_path", None, None, None
    if config.target_model_path is not None:
        if not Path(config.target_model_path).exists():
            raise FileNotFoundError(f"Target model not found: {config.target_model_path}")
        model = st.tdr.read_model(str(config.target_model_path))
        pts = np.asarray(model.points)[:, :3]
        return pts, "target_model_path", None, None, model
    if not config.build_voxel_target:
        raise ValueError("No target points supplied and build_voxel_target=False.")
    pts, pc, mesh, voxel = build_voxel_target(adata, config)
    return pts, "constructed_voxel", pc, mesh, voxel


# ---------------------------------------------------------------------------
# Interpolation methods
# ---------------------------------------------------------------------------


def run_vtk_interpolation(
    adata_dense: Any, config: InterpolationConfig,
    genes: Sequence[str], target_points: np.ndarray, **kwargs: Any,
) -> Any:
    st = _import_spateo()
    return st.tdr.vtk_interpolation(
        source_adata=adata_dense, spatial_key=config.spatial_key,
        keys=list(genes), target_points=target_points,
        n_points=config.vtk_n_points, **kwargs,
    )


def run_gp_interpolation(
    adata_dense: Any, config: InterpolationConfig,
    genes: Sequence[str], target_points: np.ndarray, **kwargs: Any,
) -> Any:
    st = _import_spateo()
    return st.tdr.gp_interpolation(
        source_adata=adata_dense, spatial_key=config.spatial_key,
        keys=list(genes), target_points=target_points,
        device=config.gp_device, **kwargs,
    )


def run_kernel_interpolation(
    adata_dense: Any, config: InterpolationConfig,
    genes: Sequence[str], target_points: np.ndarray, **kwargs: Any,
) -> Any:
    st = _import_spateo()
    return st.tdr.kernel_interpolation(
        source_adata=adata_dense, spatial_key=config.spatial_key,
        keys=list(genes), target_points=target_points, **kwargs,
    )


def run_deep_interpolation(
    adata_source: Any, config: InterpolationConfig,
    genes: Sequence[str], target_points: np.ndarray, **kwargs: Any,
) -> Any:
    st = _import_spateo()
    if hasattr(st.tdr, "deep_intepretation"):
        func = st.tdr.deep_intepretation
    elif hasattr(st.tdr, "deep_interpolation"):
        func = st.tdr.deep_interpolation
    else:
        raise AttributeError("No deep interpolation API found in st.tdr.")
    return func(
        source_adata=adata_source, spatial_key=config.spatial_key,
        keys=list(genes), target_points=target_points, **kwargs,
    )


def run_one_method(
    method: str,
    adata: Any,
    adata_dense: Any,
    config: InterpolationConfig,
    genes: Sequence[str],
    target_points: np.ndarray,
    **kwargs: Any,
) -> Any:
    method = method.lower()
    if method == "vtk":
        return run_vtk_interpolation(adata_dense, config, genes, target_points, **kwargs)
    if method == "gp":
        return run_gp_interpolation(adata_dense, config, genes, target_points, **kwargs)
    if method in {"kernel", "sparsevfc", "svfc"}:
        return run_kernel_interpolation(adata_dense, config, genes, target_points, **kwargs)
    if method == "deep":
        return run_deep_interpolation(adata, config, genes, target_points, **kwargs)
    raise ValueError(f"Unknown interpolation method: {method}")


# ---------------------------------------------------------------------------
# Gene point-cloud construction and plotting
# ---------------------------------------------------------------------------


def add_gene_labels_to_pc(pc: Any, adata: Any, genes: Sequence[str]) -> Any:
    st = _import_spateo()
    obs_index = None
    try:
        obs_index = pc.point_data["obs_index"].tolist()
    except Exception:
        pass
    for gene in genes:
        if obs_index is not None:
            try:
                x = adata[obs_index, gene].X
                exp = np.asarray(x.toarray() if hasattr(x, "toarray") else x).ravel()
            except Exception:
                exp = expression_vector(adata, gene)
        else:
            exp = expression_vector(adata, gene)
        st.tdr.add_model_labels(
            model=pc, labels=exp, key_added=gene, where="point_data", inplace=True,
        )
    return pc


def construct_gene_pc(adata: Any, spatial_key: str, genes: Sequence[str]) -> Any:
    st = _import_spateo()
    groupby = genes[0] if genes else None
    pc, _ = st.tdr.construct_pc(adata=adata.copy(), spatial_key=spatial_key, groupby=groupby)
    return add_gene_labels_to_pc(pc, adata, genes)


def save_model(model: Any, path: PathLike) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(model, "save"):
        model.save(str(path))
    else:
        st = _import_spateo()
        st.tdr.save_model(model=model, filename=str(path))
    return str(path)


def plot_gene_pc(
    pc: Any, genes: Sequence[str], filename: PathLike,
    config: InterpolationConfig, **kwargs: Any,
) -> str:
    st = _import_spateo()
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    call_kwargs: Dict[str, Any] = dict(
        model=pc, key=list(genes), colormap=config.colormap,
        opacity=config.opacity, model_style=config.model_style,
        jupyter=config.jupyter, filename=str(filename), **kwargs,
    )
    if config.camera_position is not None:
        call_kwargs["cpo"] = [config.camera_position]
    if config.window_size:
        call_kwargs["window_size"] = config.window_size
    st.pl.three_d_multi_plot(**call_kwargs)
    return str(filename)


def add_interpolated_genes_to_voxel(voxel: Any, interpolated_adata: Any, genes: Sequence[str]) -> Any:
    if voxel is None:
        return None
    result = voxel.copy() if hasattr(voxel, "copy") else voxel
    n_points = np.asarray(result.points).shape[0]
    for gene in genes:
        values = expression_vector(interpolated_adata, gene)
        if len(values) != n_points:
            LOGGER.warning("Gene %s length %d != voxel points %d; skipping.", gene, len(values), n_points)
            continue
        result.point_data[gene] = values
    return result


def plot_voxel_slices(
    voxel: Any,
    interpolated_adata: Any,
    genes: Sequence[str],
    config: InterpolationConfig,
    filename: PathLike,
    **kwargs: Any,
) -> Optional[str]:
    if voxel is None:
        LOGGER.warning("No voxel model for slice plot.")
        return None
    st = _import_spateo()
    voxel_labeled = add_interpolated_genes_to_voxel(voxel, interpolated_adata, genes)
    if voxel_labeled is None:
        return None
    LOGGER.info("Generating voxel slices along %s axis.", config.slice_axis)
    voxel_slices = st.tdr.three_d_slice(
        model=voxel_labeled, method="axis", n_slices=config.n_slices, axis=config.slice_axis,
    )
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    call_kwargs: Dict[str, Any] = dict(
        model=st.tdr.collect_models([voxel_slices]),
        key=list(genes), model_style="surface", colormap=config.colormap,
        ambient=0.5, jupyter=config.jupyter, shape=(1, len(genes)),
        text=[f"\nGene: {gene}\n" for gene in genes],
        text_kwargs={"text_size": 15, "text_font": "arial"},
        filename=str(filename), **kwargs,
    )
    if config.camera_position is not None:
        call_kwargs["cpo"] = [config.camera_position]
    st.pl.three_d_multi_plot(**call_kwargs)
    return str(filename)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def run_interpolation_pipeline(config: InterpolationConfig) -> InterpolationResult:
    configure_runtime(config)
    ensure_out_dir(config.out_dir, config.overwrite)

    adata = load_adata(config.input_path)
    validate_spatial(adata, config.spatial_key)
    adata = normalize_and_log1p(adata, config.normalize, config.log1p)
    adata = downsample_adata(
        adata, config.spatial_key, config.sample_n,
        method=config.sample_method, random_seed=config.random_seed,
    )
    genes = select_genes(adata, config.genes, config.top_n_genes)

    target_points, target_source, raw_pc, surface, voxel = get_target_points(config, adata)
    target_points = np.asarray(target_points, dtype=float)[:, :3]

    result = InterpolationResult(
        output_dir=str(config.out_dir),
        genes=list(genes),
        target_points_shape=tuple(target_points.shape),
        target_point_source=target_source,
    )

    if config.save_target_models:
        if raw_pc is not None:
            result.raw_point_cloud_model = save_model(
                raw_pc, Path(config.out_dir) / f"{config.prefix}_source_pc.vtk",
            )
        if surface is not None:
            result.surface_model = save_model(
                surface, Path(config.out_dir) / f"{config.prefix}_surface.vtk",
            )
        if voxel is not None:
            result.voxel_model = save_model(
                voxel, Path(config.out_dir) / f"{config.prefix}_voxel.vtk",
            )

    adata_dense = make_adata_dense(adata) if config.make_dense else adata

    for method in config.methods:
        method_key = method.lower()
        LOGGER.info("Running %s interpolation.", method_key)
        mr = MethodResult(method=method_key)
        try:
            interp_adata = run_one_method(
                method_key, adata, adata_dense, config, genes, target_points,
            )
            if config.save_interpolated_h5ad:
                h5ad_path = Path(config.out_dir) / f"{config.prefix}_{method_key}_interpolated.h5ad"
                interp_adata.write_h5ad(h5ad_path)
                mr.interpolated_h5ad = str(h5ad_path)
            if config.save_models or config.save_plots:
                gene_pc = construct_gene_pc(interp_adata, config.spatial_key, genes)
                if config.save_models:
                    mr.point_cloud_model = save_model(
                        gene_pc, Path(config.out_dir) / f"{config.prefix}_{method_key}_gene_pc.vtk",
                    )
                if config.save_plots:
                    mr.plot = plot_gene_pc(
                        gene_pc, genes,
                        Path(config.out_dir) / f"{config.prefix}_{method_key}_gene_expression.pdf",
                        config,
                    )
            if config.save_voxel_slices and method_key == config.slice_source_method.lower():
                mr.voxel_slice_plot = plot_voxel_slices(
                    voxel, interp_adata, genes, config,
                    Path(config.out_dir) / f"{config.prefix}_{method_key}_voxel_slices.pdf",
                )
        except Exception as exc:
            mr.error = str(exc)
            LOGGER.exception("Method %s failed.", method_key)
        result.methods[method_key] = mr

    metadata = asdict(result)
    metadata["config"] = {
        k: str(v) if isinstance(v, Path) else v for k, v in asdict(config).items()
    }
    meta_path = Path(config.out_dir) / f"{config.prefix}_interpolation_metadata.json"
    meta_path.write_text(json.dumps(metadata, indent=2, default=str, ensure_ascii=False), encoding="utf-8")
    result.metadata_path = str(meta_path)

    LOGGER.info("Interpolation complete. Output: %s", config.out_dir)
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Interpolate 3D transcriptomics signals onto reconstructed target points.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", required=True, dest="input_path")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--prefix", default="tdr_interpolation")
    p.add_argument("--spatial-key", default="3d_align_spatial")
    p.add_argument("--groupby", default=None)
    p.add_argument("--genes", nargs="+", default=[])
    p.add_argument("--top-n-genes", type=int, default=0)
    p.add_argument("--normalize", action="store_true")
    p.add_argument("--log1p", action="store_true")
    p.add_argument("--no-dense", action="store_true")
    p.add_argument("--sample-n", type=int, default=None)
    p.add_argument("--sample-method", default="random")
    p.add_argument("--random-seed", type=int, default=0)
    p.add_argument("--methods", nargs="+", default=["vtk", "gp", "kernel"],
                   choices=["vtk", "gp", "kernel", "sparsevfc", "svfc", "deep"])
    p.add_argument("--target-points", default=None, dest="target_points_path")
    p.add_argument("--target-model", default=None, dest="target_model_path")
    p.add_argument("--target-point-columns", nargs=3, default=("x", "y", "z"))
    p.add_argument("--no-build-voxel-target", action="store_true")
    p.add_argument("--surface-alpha", type=float, default=0.6)
    p.add_argument("--surface-method", default="marching_cube")
    p.add_argument("--surface-mc-scale-factor", type=float, default=1.0)
    p.add_argument("--surface-smooth", type=int, default=5000)
    p.add_argument("--surface-scale-factor", type=float, default=1.0)
    p.add_argument("--voxel-smooth", type=int, default=50)
    p.add_argument("--vtk-n-points", type=int, default=5)
    p.add_argument("--gp-device", default="cpu")
    p.add_argument("--no-save-h5ad", action="store_true")
    p.add_argument("--no-save-models", action="store_true")
    p.add_argument("--no-save-target-models", action="store_true")
    p.add_argument("--save-plots", action="store_true")
    p.add_argument("--save-voxel-slices", action="store_true")
    p.add_argument("--slice-axis", default="x", choices=["x", "y", "z"])
    p.add_argument("--n-slices", type=int, default=15)
    p.add_argument("--slice-source-method", default="gp")
    p.add_argument("--colormap", default="hot_r")
    p.add_argument("--opacity", type=float, default=0.5)
    p.add_argument("--model-style", default="points")
    p.add_argument("--jupyter", default="static")
    p.add_argument("--window-size", nargs=2, type=int, default=[500, 500])
    p.add_argument("--camera-position", default=None)
    p.add_argument("--xvfb", action="store_true", dest="start_xvfb")
    p.add_argument("--off-screen", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    configure_logging(verbose=args.verbose)
    try:
        config = InterpolationConfig(
            input_path=args.input_path,
            out_dir=args.out_dir,
            prefix=args.prefix,
            spatial_key=args.spatial_key,
            groupby=args.groupby,
            genes=args.genes,
            top_n_genes=args.top_n_genes,
            normalize=args.normalize,
            log1p=args.log1p,
            make_dense=not args.no_dense,
            sample_n=args.sample_n,
            sample_method=args.sample_method,
            random_seed=args.random_seed,
            methods=args.methods,
            target_points_path=args.target_points_path,
            target_model_path=args.target_model_path,
            target_point_columns=tuple(args.target_point_columns),
            build_voxel_target=not args.no_build_voxel_target,
            surface_alpha=args.surface_alpha,
            surface_method=args.surface_method,
            surface_mc_scale_factor=args.surface_mc_scale_factor,
            surface_smooth=args.surface_smooth,
            surface_scale_factor=args.surface_scale_factor,
            voxel_smooth=args.voxel_smooth,
            vtk_n_points=args.vtk_n_points,
            gp_device=args.gp_device,
            save_interpolated_h5ad=not args.no_save_h5ad,
            save_models=not args.no_save_models,
            save_target_models=not args.no_save_target_models,
            save_plots=args.save_plots,
            save_voxel_slices=args.save_voxel_slices,
            slice_axis=args.slice_axis,
            n_slices=args.n_slices,
            slice_source_method=args.slice_source_method,
            colormap=args.colormap,
            opacity=args.opacity,
            model_style=args.model_style,
            jupyter=args.jupyter,
            window_size=(args.window_size[0], args.window_size[1]),
            camera_position=parse_camera_position(args.camera_position),
            start_xvfb=args.start_xvfb,
            off_screen=args.off_screen,
            overwrite=args.overwrite,
        )
        result = run_interpolation_pipeline(config)
    except Exception as exc:
        LOGGER.exception("Interpolation pipeline failed: %s", exc)
        return 1
    print(json.dumps(asdict(result), indent=2, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "InterpolationConfig",
    "InterpolationResult",
    "MethodResult",
    "build_arg_parser",
    "configure_logging",
    "main",
    "run_deep_interpolation",
    "run_gp_interpolation",
    "run_interpolation_pipeline",
    "run_kernel_interpolation",
    "run_one_method",
    "run_vtk_interpolation",
    "validate_spatial",
]
