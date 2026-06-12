"""Stage 3 morphogenesis: cell mapping, vector field, and trajectory analysis."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

PathLike = Union[str, Path]

LOGGER = logging.getLogger("spateo_skills.morphogenesis.vectorfield")

DEFAULT_CPO = [
    (80.14473968731744, 416.994398326684, 633.0950009292029),
    (2.2346878051757812, -11.286853790283203, 1.3404006958007812),
    (-0.29255021850853963, -0.774368934684527, 0.5610411060217927),
]


@dataclass
class VectorFieldConfig:
    stage1_path: str = ""
    stage2_path: str = ""
    output_path: Optional[str] = None
    preprocess: bool = False
    compute_trajectory: bool = True
    device: str = "0"
    counts_layer: str = "counts_X"
    spatial_key: str = "3d_align_spatial"
    layer: str = "log1p_X"
    vf_key: str = "VecFld_morpho"
    fate_key: str = "fate_morpho"
    mapping_key_added: str = "cells_mapping"
    numItermaxEmd: int = 500000
    t_end: int = 20000
    interpolation_num: int = 50
    cores: int = 20
    stage1_mc_scale_factor: float = 1.22
    stage2_mc_scale_factor: float = 1.12
    scale_factor: float = 1.08
    save_prefix: str = "amn"


@dataclass
class VectorFieldResult:
    vf_key: str = ""
    saved_model_paths: List[str] = field(default_factory=list)
    saved_adata_path: str = ""
    saved_trajectory_path: str = ""
    pot_patched: bool = False
    notes: List[str] = field(default_factory=list)


def patch_pot_cg_if_needed() -> bool:
    """Apply POT compatibility patch; return True if applied."""
    try:
        import ot
    except Exception:
        return False
    if not hasattr(ot.gromov, "cg") and hasattr(ot.optim, "cg"):
        ot.gromov.cg = ot.optim.cg
        return True
    return False


def _nonzero_gene_mask(layer: Any) -> np.ndarray:
    sums = np.asarray(layer.sum(axis=0)).ravel()
    return sums != 0


def load_stage_adata(
    path: PathLike,
    *,
    counts_layer: str = "counts_X",
    reset_pp: bool = True,
    set_X_to_counts: bool = True,
    filter_zero_genes: bool = True,
) -> Any:
    """Load and minimally normalize one stage AnnData."""
    import spateo as st

    adata = st.read_h5ad(str(Path(path)))
    if filter_zero_genes and counts_layer in adata.layers:
        adata = adata[:, _nonzero_gene_mask(adata.layers[counts_layer])].copy()
    if reset_pp:
        adata.uns["pp"] = {}
    if set_X_to_counts and counts_layer in adata.layers:
        adata.X = adata.layers[counts_layer].copy()
    return adata


def load_stage_pair(
    stage1_path: PathLike,
    stage2_path: PathLike,
    **kwargs: Any,
) -> Tuple[Any, Any]:
    """Load two stages for morphogenesis analysis."""
    return load_stage_adata(stage1_path, **kwargs), load_stage_adata(stage2_path, **kwargs)


def preprocess_stage_pair(
    stage1_adata: Any,
    stage2_adata: Any,
    *,
    recipe: str = "generic",
    spatial_key: str = "spatial",
    counts_layer: str = "counts",
    min_genes: int = 10,
    min_cells: int = 3,
    feature_method: str = "hvg",
    n_top_genes: int = 2000,
    **kwargs: Any,
) -> Tuple[Any, Any]:
    """Run preprocess_spatial on both stages."""
    from spateo.preprocessing.protocol_pipeline import preprocess_spatial

    for adata in (stage1_adata, stage2_adata):
        preprocess_spatial(
            adata,
            recipe=recipe,
            spatial_key=spatial_key,
            counts_layer=counts_layer,
            min_genes=min_genes,
            min_cells=min_cells,
            feature_method=feature_method,
            n_top_genes=n_top_genes,
            **kwargs,
        )
    return stage1_adata, stage2_adata


def compute_cell_mapping(
    stage1_adata: Any,
    stage2_adata: Any,
    *,
    layer: str = "log1p_X",
    spatial_key: str = "3d_align_spatial",
    key_added: str = "cells_mapping",
    alpha: float = 0.0001,
    numItermaxEmd: int = 500000,
    device: Union[str, int] = "0",
    inplace: bool = True,
    **kwargs: Any,
) -> Any:
    """Map cells from stage A to stage B via spateo.tdr.cell_directions."""
    import spateo as st

    return st.tdr.cell_directions(
        adataA=stage1_adata,
        adataB=stage2_adata,
        layer=layer,
        numItermaxEmd=numItermaxEmd,
        spatial_key=spatial_key,
        key_added=key_added,
        alpha=alpha,
        device=device,
        inplace=inplace,
        **kwargs,
    )


def construct_stage_pc_and_mesh(
    adata: Any,
    *,
    spatial_key: str = "3d_align_spatial",
    groupby: str = "Annotation_2_tissue",
    key_added: str = "tissue",
    colormap: Optional[Dict[str, str]] = None,
    alpha: float = 0.6,
    mc_scale_factor: float = 1.2,
    smooth: int = 8000,
    scale_factor: float = 1.08,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Construct point cloud and surface mesh for one tissue/stage."""
    import spateo as st

    if colormap is None:
        colormap = {"amnioserosa": "#DC143C"}
    pc, _ = st.tdr.construct_pc(
        adata=adata,
        spatial_key=spatial_key,
        groupby=groupby,
        key_added=key_added,
        colormap=colormap,
    )
    mesh, _, _ = st.tdr.construct_surface(
        pc=pc,
        key_added=key_added,
        alpha=alpha,
        cs_method="marching_cube",
        cs_args={"mc_scale_factor": mc_scale_factor},
        smooth=smooth,
        scale_factor=scale_factor,
    )
    return {"pc": pc, "mesh": mesh}


def construct_stage_models(
    stage1_adata: Any,
    stage2_adata: Any,
    *,
    stage1_mc_scale_factor: float = 1.22,
    stage2_mc_scale_factor: float = 1.12,
    scale_factor: float = 1.08,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Construct pc/mesh models for both stages."""
    stage1 = construct_stage_pc_and_mesh(
        stage1_adata,
        mc_scale_factor=stage1_mc_scale_factor,
        scale_factor=scale_factor,
        **kwargs,
    )
    stage2 = construct_stage_pc_and_mesh(
        stage2_adata,
        mc_scale_factor=stage2_mc_scale_factor,
        scale_factor=scale_factor,
        **kwargs,
    )
    return {
        "stage1_pc": stage1["pc"],
        "stage1_mesh": stage1["mesh"],
        "stage2_pc": stage2["pc"],
        "stage2_mesh": stage2["mesh"],
    }


def save_models_and_adata(
    *,
    outdir: PathLike,
    stage1_pc: Any,
    stage1_mesh: Any,
    stage2_pc: Optional[Any] = None,
    stage2_mesh: Optional[Any] = None,
    stage1_adata: Optional[Any] = None,
    prefix: str = "amn",
    **kwargs: Any,
) -> Dict[str, str]:
    """Save models and stage1 AnnData to disk; return saved paths."""
    import spateo as st

    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    saved: Dict[str, str] = {}
    model_specs = {
        "stage1_aligned_pc": stage1_pc,
        "stage1_aligned_mesh": stage1_mesh,
        "stage2_aligned_pc": stage2_pc,
        "stage2_aligned_mesh": stage2_mesh,
    }
    for name, model in model_specs.items():
        if model is None:
            continue
        path = outdir_path / f"{prefix}_{name}.vtk"
        st.tdr.save_model(model=model, filename=str(path))
        saved[name] = str(path)
    if stage1_adata is not None:
        path = outdir_path / f"{prefix}_stage1_unmorphed.h5ad"
        stage1_adata.write(str(path))
        saved["stage1_adata"] = str(path)
    return saved


def construct_mapping_lines(
    stage1_adata: Any,
    *,
    spatial_key: str = "3d_align_spatial",
    mapping_key: str = "X_cells_mapping",
    z_offset: float = -600.0,
    key_added: str = "check_align",
    label: str = "align_lines",
    color: str = "gainsboro",
    **kwargs: Any,
) -> Any:
    """Construct lines visualizing cell mappings between stages."""
    import spateo as st

    return st.tdr.construct_align_lines(
        model1_points=stage1_adata.obsm[spatial_key].copy(),
        model2_points=stage1_adata.obsm[mapping_key].copy() + np.asarray([0, 0, z_offset]),
        key_added=key_added,
        label=label,
        color=color,
    )


def compute_morphogenesis_vectorfield(
    stage1_adata: Any,
    stage1_mesh: Any,
    *,
    spatial_key: str = "3d_align_spatial",
    V_key: str = "V_cells_mapping",
    key_added: str = "VecFld_morpho",
    inplace: bool = True,
    **kwargs: Any,
) -> Any:
    """Compute developmental vector field via morphofield_sparsevfc."""
    import spateo as st

    return st.tdr.morphofield_sparsevfc(
        adata=stage1_adata,
        spatial_key=spatial_key,
        V_key=V_key,
        key_added=key_added,
        NX=np.asarray(stage1_mesh.points),
        inplace=inplace,
        **kwargs,
    )


def attach_vectorfield_to_models(
    stage1_adata: Any,
    stage1_pc: Any,
    stage1_mesh: Any,
    *,
    vf_key: str = "VecFld_morpho",
    model_vector_key: str = "vectors",
) -> Tuple[Any, Any]:
    """Attach computed vector-field values to point cloud and mesh models."""
    stage1_pc.point_data[model_vector_key] = stage1_adata.uns[vf_key]["V"]
    stage1_mesh.point_data[model_vector_key] = stage1_adata.uns[vf_key]["grid_V"]
    return stage1_pc, stage1_mesh


def construct_vector_arrows(
    model: Any,
    *,
    vf_key: str = "vectors",
    factor: float = 20000.0,
    key_added: str = "v_arrows",
    label_component: int = 2,
    color: str = "rainbow",
    **kwargs: Any,
) -> Any:
    """Construct arrow glyphs for a vector field attached to a model."""
    import spateo as st

    return st.tdr.construct_field(
        model=model,
        vf_key=vf_key,
        arrows_scale_key=vf_key,
        factor=factor,
        key_added=key_added,
        label=model.point_data[vf_key][:, label_component].flatten(),
        color=color,
        **kwargs,
    )


def compute_developmental_trajectory(
    stage1_adata: Any,
    *,
    layer: str = "log1p_X",
    vf_key: str = "VecFld_morpho",
    key_added: str = "fate_morpho",
    t_end: int = 20000,
    interpolation_num: int = 50,
    cores: int = 20,
    sampling_method: str = "trn",
    trajectory_color: str = "rainbow",
    tip_color: str = "rainbow",
    **kwargs: Any,
) -> Any:
    """Predict developmental trajectories via morphopath + construct_trajectory."""
    import spateo as st

    st.tdr.morphopath(
        adata=stage1_adata,
        layer=layer,
        vf_key=vf_key,
        key_added=key_added,
        t_end=t_end,
        interpolation_num=interpolation_num,
        cores=cores,
    )
    trajectory_model, _ = st.tdr.construct_trajectory(
        adata=stage1_adata,
        fate_key=key_added,
        sampling_method=sampling_method,
        label=stage1_adata.uns[vf_key]["V"][:, 2].flatten(),
        trajectory_color=trajectory_color,
        tip_color=tip_color,
        **kwargs,
    )
    return trajectory_model


def plot_models(model: Any, *, key: Any = "tissue", cpo: Any = DEFAULT_CPO, **kwargs: Any) -> Any:
    """Generic wrapper around spateo.pl.three_d_plot."""
    import spateo as st

    defaults = dict(jupyter="static", show_axes=True, cpo=cpo)
    defaults.update(kwargs)
    return st.pl.three_d_plot(model=model, key=key, **defaults)


def run_vectorfield_pipeline(config: VectorFieldConfig, **kwargs: Any) -> VectorFieldResult:
    """End-to-end morphogenesis mapping, vector-field, and trajectory pipeline."""
    result = VectorFieldResult(vf_key=config.vf_key)
    result.pot_patched = patch_pot_cg_if_needed()

    stage1_adata, stage2_adata = load_stage_pair(
        config.stage1_path,
        config.stage2_path,
        counts_layer=config.counts_layer,
    )
    result.notes.append(f"Loaded stage1={config.stage1_path}, stage2={config.stage2_path}")

    if config.preprocess:
        preprocess_stage_pair(stage1_adata, stage2_adata, **kwargs)
        result.notes.append("Preprocessed both stages")

    compute_cell_mapping(
        stage1_adata,
        stage2_adata,
        layer=config.layer,
        spatial_key=config.spatial_key,
        key_added=config.mapping_key_added,
        numItermaxEmd=config.numItermaxEmd,
        device=config.device,
        **kwargs,
    )
    result.notes.append("Computed cell mapping")

    models = construct_stage_models(
        stage1_adata,
        stage2_adata,
        stage1_mc_scale_factor=config.stage1_mc_scale_factor,
        stage2_mc_scale_factor=config.stage2_mc_scale_factor,
        scale_factor=config.scale_factor,
        **kwargs,
    )
    result.notes.append("Constructed stage models (pc + mesh)")

    compute_morphogenesis_vectorfield(
        stage1_adata,
        models["stage1_mesh"],
        spatial_key=config.spatial_key,
        key_added=config.vf_key,
        **kwargs,
    )
    result.notes.append(f"Computed morphogenesis vectorfield (key={config.vf_key})")

    attach_vectorfield_to_models(
        stage1_adata,
        models["stage1_pc"],
        models["stage1_mesh"],
        vf_key=config.vf_key,
    )

    if config.compute_trajectory:
        compute_developmental_trajectory(
            stage1_adata,
            layer=config.layer,
            vf_key=config.vf_key,
            key_added=config.fate_key,
            t_end=config.t_end,
            interpolation_num=config.interpolation_num,
            cores=config.cores,
            **kwargs,
        )
        result.notes.append(f"Computed developmental trajectory (key={config.fate_key})")

    if config.output_path:
        saved = save_models_and_adata(
            outdir=config.output_path,
            stage1_pc=models["stage1_pc"],
            stage1_mesh=models["stage1_mesh"],
            stage2_pc=models["stage2_pc"],
            stage2_mesh=models["stage2_mesh"],
            stage1_adata=stage1_adata,
            prefix=config.save_prefix,
        )
        result.saved_model_paths = [
            v for k, v in saved.items() if k != "stage1_adata"
        ]
        result.saved_adata_path = saved.get("stage1_adata", "")
        result.notes.append(f"Saved models to {config.output_path}")

    return result


def configure_logging(verbose: bool = False) -> None:
    """Configure root logger for this module."""
    from ..shared.data_helpers import configure_logging as _configure

    _configure(verbose)


def build_arg_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser for vectorfield pipeline."""
    parser = argparse.ArgumentParser(
        description="Morphogenesis vectorfield pipeline: cell mapping, vector field, trajectory.",
    )
    parser.add_argument("--stage1", required=True, help="Path to stage 1 h5ad")
    parser.add_argument("--stage2", required=True, help="Path to stage 2 h5ad")
    parser.add_argument("--output", default=None, help="Output directory for models")
    parser.add_argument("--preprocess", action="store_true", help="Run preprocessing on both stages")
    parser.add_argument(
        "--no-trajectory",
        action="store_false",
        dest="compute_trajectory",
        help="Skip trajectory computation",
    )
    parser.add_argument("--device", default="0", help="CUDA device index or 'cpu'")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    parser.add_argument("--prefix", default="amn", help="Save file prefix")
    parser.set_defaults(compute_trajectory=True)
    return parser


def main() -> int:
    """CLI entry point."""
    parser = build_arg_parser()
    args = parser.parse_args()

    configure_logging(args.verbose)

    config = VectorFieldConfig(
        stage1_path=args.stage1,
        stage2_path=args.stage2,
        output_path=args.output,
        preprocess=args.preprocess,
        compute_trajectory=args.compute_trajectory,
        device=args.device,
        save_prefix=args.prefix,
    )

    result = run_vectorfield_pipeline(config)
    print(json.dumps(asdict(result), indent=2, ensure_ascii=False))
    return 0


__all__ = [
    "DEFAULT_CPO",
    "VectorFieldConfig",
    "VectorFieldResult",
    "patch_pot_cg_if_needed",
    "load_stage_adata",
    "load_stage_pair",
    "preprocess_stage_pair",
    "compute_cell_mapping",
    "construct_stage_pc_and_mesh",
    "construct_stage_models",
    "save_models_and_adata",
    "construct_mapping_lines",
    "compute_morphogenesis_vectorfield",
    "attach_vectorfield_to_models",
    "construct_vector_arrows",
    "compute_developmental_trajectory",
    "plot_models",
    "run_vectorfield_pipeline",
    "configure_logging",
    "build_arg_parser",
    "main",
]
