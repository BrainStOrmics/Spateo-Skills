"""Stage 3 morphogenesis: feature computation and GLM DEG analysis.

Extracts dynamic features from vector fields (velocity, acceleration,
curvature, curl, torsion, jacobian) and identifies differentially
expressed genes via GLM.

Usage:
    python -m skills.04_morphogenesis.feature 
        --input ./data/vectorfield.h5ad --out-dir ./output
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

PathLike = Union[str, Path]

LOGGER = logging.getLogger("spateo_skills.morphogenesis.feature")

DEFAULT_CPO = [
    (80.14473968731744, 416.994398326684, 633.0950009292029),
    (2.2346878051757812, -11.286853790283203, 1.3404006958007812),
    (-0.29255021850853963, -0.774368934684527, 0.5610411060217927),
]

FEATURE_FUNCTIONS = {
    "velocity": "morphofield_velocity",
    "acceleration": "morphofield_acceleration",
    "curvature": "morphofield_curvature",
    "curl": "morphofield_curl",
    "torsion": "morphofield_torsion",
    "jacobian": "morphofield_jacobian",
}

PLOT_FUNCTIONS = {
    "acceleration": "acceleration",
    "curvature": "curvature",
    "curl": "curl",
    "torsion": "torsion",
    "jacobian": "jacobian",
}


@dataclass
class FeatureConfig:
    adata_path: str = ""
    pc_model_path: Optional[str] = None
    mesh_model_path: Optional[str] = None
    output_path: Optional[str] = None
    features: List[str] = field(
        default_factory=lambda: list(FEATURE_FUNCTIONS.keys()),
    )
    compute_trajectory: bool = True
    run_glm: bool = True
    counts_layer: str = "counts_X"
    vf_key: str = "VecFld_morpho"
    fate_key: str = "fate_morpho"


@dataclass
class FeatureResult:
    feature_keys: Dict[str, str] = field(default_factory=dict)
    glm_keys: Dict[str, str] = field(default_factory=dict)
    saved_table_paths: Dict[str, str] = field(default_factory=dict)
    trajectory_model_key: str = ""
    notes: List[str] = field(default_factory=list)


def _nonzero_gene_mask(layer: Any) -> np.ndarray:
    sums = np.asarray(layer.sum(axis=0)).ravel()
    return sums != 0


def load_morphogenesis_inputs(
    adata_path: PathLike,
    pc_model_path: Optional[PathLike] = None,
    mesh_model_path: Optional[PathLike] = None,
    *,
    counts_layer: str = "counts_X",
    set_X_to_counts: bool = True,
    filter_zero_genes: bool = True,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Load AnnData and optional pc/mesh models from vectorfield stage."""
    import spateo as st

    adata = st.read_h5ad(str(Path(adata_path)))
    if filter_zero_genes and counts_layer in adata.layers:
        adata = adata[:, _nonzero_gene_mask(adata.layers[counts_layer])].copy()
    if set_X_to_counts and counts_layer in adata.layers:
        adata.X = adata.layers[counts_layer].copy()

    pc = st.tdr.read_model(str(Path(pc_model_path))) if pc_model_path is not None else None
    mesh = st.tdr.read_model(str(Path(mesh_model_path))) if mesh_model_path is not None else None
    return {"adata": adata, "pc": pc, "mesh": mesh}


def compute_or_load_trajectory(
    adata: Any,
    *,
    layer: str = "log1p_X",
    vf_key: str = "VecFld_morpho",
    fate_key: str = "fate_morpho",
    t_end: int = 20000,
    interpolation_num: int = 50,
    cores: int = 20,
    key_added: str = "obs_index",
    label_from_obs_index: bool = True,
    **kwargs: Any,
) -> Any:
    """Compute morphopath and construct a trajectory model."""
    import spateo as st

    st.tdr.morphopath(
        adata=adata,
        layer=layer,
        vf_key=vf_key,
        key_added=fate_key,
        t_end=t_end,
        interpolation_num=interpolation_num,
        cores=cores,
    )
    label = np.asarray(adata.obs.index) if label_from_obs_index else None
    trajectory_model, _ = st.tdr.construct_trajectory(
        adata=adata,
        fate_key=fate_key,
        key_added=key_added,
        label=label,
        **kwargs,
    )
    return trajectory_model


def compute_morphogenesis_feature(
    adata: Any,
    feature: str,
    *,
    vf_key: str = "VecFld_morpho",
    key_added: Optional[str] = None,
    **kwargs: Any,
) -> str:
    """Compute one morphogenesis feature; return the key where it was stored."""
    import spateo as st

    feature = feature.lower()
    if feature not in FEATURE_FUNCTIONS:
        raise ValueError(f"Unsupported feature {feature!r}; supported: {sorted(FEATURE_FUNCTIONS)}")
    key = key_added or feature
    func = getattr(st.tdr, FEATURE_FUNCTIONS[feature])
    func(adata=adata, vf_key=vf_key, key_added=key, **kwargs)
    return key


def run_glm_for_feature(
    adata: Any,
    feature_key: str,
    *,
    glm_key: Optional[str] = None,
    layer: Optional[str] = None,
    qval_threshold: float = 0.25,
    llf_threshold: Optional[float] = None,
    df: int = 3,
    **kwargs: Any,
) -> tuple:
    """Run GLM differential gene analysis against a morphogenesis feature."""
    import spateo as st

    glm_key = glm_key or f"glm_degs_{feature_key}"
    st.tl.glm_degs(
        adata=adata,
        layer=layer,
        fullModelFormulaStr=f"~cr({feature_key}, df={df})",
        key_added=glm_key,
        qval_threshold=qval_threshold,
        llf_threshold=llf_threshold,
        **kwargs,
    )
    return adata.uns[glm_key]["glm_result"], glm_key


def compute_feature_with_glm(
    adata: Any,
    feature: str,
    *,
    vf_key: str = "VecFld_morpho",
    qval_threshold: Optional[float] = None,
    llf_threshold: Optional[float] = None,
    run_glm: bool = True,
    glm_dict: Optional[Dict[str, Any]] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Compute a feature and optionally its GLM result."""
    feature = feature.lower()
    feature_key = compute_morphogenesis_feature(adata, feature, vf_key=vf_key)
    out: Dict[str, Any] = {"feature": feature, "feature_key": feature_key}
    if run_glm and feature not in ("velocity", "jacobian"):
        default_q = {
            "acceleration": 0.01,
            "curvature": 0.25,
            "curl": 0.25,
            "torsion": 0.01,
        }.get(feature, 0.25)
        default_llf: Optional[float] = {"curl": -50.0}.get(feature, None)
        glm_data, glm_key = run_glm_for_feature(
            adata,
            feature_key,
            qval_threshold=default_q if qval_threshold is None else qval_threshold,
            llf_threshold=default_llf if llf_threshold is None else llf_threshold,
            **kwargs,
        )
        out.update({"glm_key": glm_key, "glm_result": glm_data})
        if glm_dict is not None:
            glm_dict[feature] = glm_data
    return out


def plot_glm_fit(
    adata: Any,
    glm_data: Any,
    *,
    feature_key: str,
    glm_key: str,
    ncols: int = 2,
    save_show_or_return: str = "return",
    **kwargs: Any,
) -> Any:
    """Plot GLM fitted genes for one feature."""
    import spateo as st

    return st.pl.glm_fit(
        adata=adata,
        genes=glm_data.index.tolist(),
        ncols=ncols,
        feature_x=feature_key,
        feature_y="expression",
        glm_key=glm_key,
        save_show_or_return=save_show_or_return,
        **kwargs,
    )


def plot_morphogenesis_feature(
    adata: Any,
    feature: str,
    model: Any,
    *,
    key: Optional[str] = None,
    cpo: Any = DEFAULT_CPO,
    filename: Optional[PathLike] = None,
    **kwargs: Any,
) -> Any:
    """Plot acceleration/curvature/curl/torsion/jacobian on a model."""
    import spateo as st

    feature = feature.lower()
    if feature not in PLOT_FUNCTIONS:
        raise ValueError(f"No plot wrapper for feature {feature!r}; supported: {sorted(PLOT_FUNCTIONS)}")
    key = key or feature
    plot_func = getattr(st.pl, PLOT_FUNCTIONS[feature])
    defaults = dict(
        adata=adata,
        model=model,
        colormap="default_cmap",
        jupyter="static",
        model_style=["points", "wireframe"],
        model_size=[5, 2],
        cpo=cpo,
        window_size=(1024, 1024),
        background="white",
    )
    defaults.update(kwargs)
    if filename is not None:
        defaults["filename"] = str(filename)

    if feature == "jacobian":
        defaults["jacobian_key"] = key
    else:
        defaults[f"{feature}_key"] = key
    return plot_func(**defaults)


def run_feature_pipeline(config: FeatureConfig, **kwargs: Any) -> FeatureResult:
    """End-to-end feature computation pipeline."""
    result = FeatureResult()

    inputs = load_morphogenesis_inputs(
        config.adata_path,
        config.pc_model_path,
        config.mesh_model_path,
        counts_layer=config.counts_layer,
    )
    adata = inputs["adata"]
    result.notes.append(f"Loaded adata from {config.adata_path}")

    if config.compute_trajectory:
        compute_or_load_trajectory(
            adata,
            vf_key=config.vf_key,
            fate_key=config.fate_key,
            **kwargs,
        )
        result.trajectory_model_key = config.fate_key
        result.notes.append(f"Computed trajectory (key={config.fate_key})")

    glm_dict: Dict[str, Any] = {}
    for feature in config.features:
        feature_out = compute_feature_with_glm(
            adata,
            feature,
            vf_key=config.vf_key,
            run_glm=config.run_glm,
            glm_dict=glm_dict,
            **kwargs,
        )
        result.feature_keys[feature] = feature_out["feature_key"]
        if "glm_key" in feature_out:
            result.glm_keys[feature] = feature_out["glm_key"]
        result.notes.append(f"Computed feature: {feature}")

    if config.output_path:
        outdir = Path(config.output_path)
        outdir.mkdir(parents=True, exist_ok=True)
        for feature, glm_data in glm_dict.items():
            path = outdir / f"glm_degs_{feature}.csv"
            glm_data.to_csv(path)
            result.saved_table_paths[feature] = str(path)
        result.notes.append(f"Saved GLM tables to {config.output_path}")

    return result


def configure_logging(verbose: bool = False) -> None:
    """Configure root logger for this module."""
    from skills.shared.data_helpers import configure_logging as _configure

    _configure(verbose)


def build_arg_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser for feature pipeline."""
    parser = argparse.ArgumentParser(
        description="Morphogenesis feature pipeline: velocity/acceleration/curvature/curl/torsion/jacobian + GLM DEG.",
    )
    parser.add_argument("--adata", required=True, help="Path to AnnData h5ad (from vectorfield stage)")
    parser.add_argument("--pc-model", default=None, help="Path to point cloud VTK model")
    parser.add_argument("--mesh-model", default=None, help="Path to mesh VTK model")
    parser.add_argument("--output", default=None, help="Output directory for GLM CSV tables")
    parser.add_argument(
        "--features",
        nargs="+",
        default=list(FEATURE_FUNCTIONS.keys()),
        help="Features to compute",
    )
    parser.add_argument(
        "--no-trajectory",
        action="store_false",
        dest="compute_trajectory",
        help="Skip trajectory computation",
    )
    parser.add_argument(
        "--no-glm",
        action="store_false",
        dest="run_glm",
        help="Skip GLM DEG analysis",
    )
    parser.add_argument("--vf-key", default="VecFld_morpho", help="Vector field key in adata.uns")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    parser.set_defaults(compute_trajectory=True, run_glm=True)
    return parser


def main() -> int:
    """CLI entry point."""
    parser = build_arg_parser()
    args = parser.parse_args()

    configure_logging(args.verbose)

    config = FeatureConfig(
        adata_path=args.adata,
        pc_model_path=args.pc_model,
        mesh_model_path=args.mesh_model,
        output_path=args.output,
        features=args.features,
        compute_trajectory=args.compute_trajectory,
        run_glm=args.run_glm,
        vf_key=args.vf_key,
    )

    result = run_feature_pipeline(config)
    print(json.dumps(asdict(result), indent=2, ensure_ascii=False))
    return 0


__all__ = [
    "DEFAULT_CPO",
    "FEATURE_FUNCTIONS",
    "PLOT_FUNCTIONS",
    "FeatureConfig",
    "FeatureResult",
    "load_morphogenesis_inputs",
    "compute_or_load_trajectory",
    "compute_morphogenesis_feature",
    "run_glm_for_feature",
    "compute_feature_with_glm",
    "plot_glm_fit",
    "plot_morphogenesis_feature",
    "run_feature_pipeline",
    "configure_logging",
    "build_arg_parser",
    "main",
]
