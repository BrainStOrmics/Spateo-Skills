#!/usr/bin/env python3
"""run_all — Skill: Full Spateo pipeline orchestrator.

Runs Stage 0 (IO) → Stage 1 (Alignment) → Stage 2 (Reconstruction) → Stage 3 (Morphogenesis)
in sequence, with per-stage skip flags and shared configuration.

Usage:
    python -m skills.00_pipeline.run_all \
        --data-path ./data/xenium_outs/ \
        --platform xenium \
        --stage2-path ./data/stage2.h5ad \
        --out-dir ./output/
"""
from __future__ import annotations

import argparse
import importlib
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional, Sequence

from ..shared.data_helpers import (
    configure_logging,
    infer_device,
    load_h5ad,
)

LOGGER = logging.getLogger("pipeline.run_all")


@dataclass
class PipelineConfig:
    data_path: Optional[str] = None
    platform: Optional[str] = None
    auto_detect: bool = False
    slices_dir: Optional[str] = None
    stage2_path: Optional[str] = None
    out_dir: str = "./output/"
    skip_alignment: bool = False
    skip_reconstruction: bool = False
    skip_morphogenesis: bool = False
    model_type: str = "surface"
    interpolation_method: str = "vtk"
    groupby: str = "tissue_type"
    device: Optional[str] = None
    verbose: bool = False
    stage0_kwargs: dict = field(default_factory=dict)
    stage1_kwargs: dict = field(default_factory=dict)
    stage2_kwargs: dict = field(default_factory=dict)
    stage3_kwargs: dict = field(default_factory=dict)


@dataclass
class PipelineResult:
    stage0_output: str = ""
    stage1_output: str = ""
    stage2_output: str = ""
    stage3_output: str = ""
    notes: list[str] = field(default_factory=list)


def run_stage0(config: PipelineConfig) -> Any:
    auto_read = importlib.import_module("skills.01_data_io.auto_read")
    manual_read = importlib.import_module("skills.01_data_io.manual_read")

    if config.auto_detect:
        adata, match = auto_read.read_auto_spatial_data(config.data_path, return_match=True)
        LOGGER.info("Auto-detected platform: %s", match)
        return adata
    elif config.platform:
        return manual_read.read_by_platform(config.platform, config.data_path, **config.stage0_kwargs)
    else:
        raise ValueError("Set --platform or --auto-detect for Stage 0")


def run_stage1(config: PipelineConfig, stage0_adata: Any) -> Any:
    two_slice = importlib.import_module("skills.02_slice_alignment.two_slice")

    if config.slices_dir:
        slices = sorted(Path(config.slices_dir).glob("*.h5ad"))
        if len(slices) < 2:
            raise ValueError(f"Need >= 2 .h5ad files in {config.slices_dir}")
        cfg = two_slice.TwoSliceConfig(
            slice1_path=str(slices[0]),
            slice2_path=str(slices[1]),
            **config.stage1_kwargs,
        )
        return two_slice.run_two_slice_pipeline(cfg)
    else:
        raise ValueError("Set --slices-dir or provide slice paths for Stage 1")


def run_stage2(config: PipelineConfig, aligned_adata: Any) -> str:
    reconstruction = importlib.import_module("skills.03_3d_reconstruction.reconstruction")

    tmp_path = Path(config.out_dir) / "_stage2_aligned.h5ad"
    tmp_path.parent.mkdir(parents=True, exist_ok=True)
    aligned_adata.write_h5ad(tmp_path)

    cfg = reconstruction.ReconstructionConfig(
        input_path=str(tmp_path),
        model_type=config.model_type,
        groupby=config.groupby,
        out_dir=config.out_dir,
        **config.stage2_kwargs,
    )
    return reconstruction.run_reconstruction_pipeline(cfg).outdir


def run_stage3(config: PipelineConfig, stage2_adata: Any) -> str:
    if not config.stage2_path:
        raise ValueError("Set --stage2-path for Stage 3 (morphogenesis)")
    vectorfield = importlib.import_module("skills.04_morphogenesis.vectorfield")

    cfg = vectorfield.VectorFieldConfig(
        stage1_path=str(Path(config.out_dir) / "_stage2_aligned.h5ad"),
        stage2_path=config.stage2_path,
        groupby=config.groupby,
        **config.stage3_kwargs,
    )
    return vectorfield.run_vectorfield_pipeline(cfg).outdir


def run_full_pipeline(config: PipelineConfig) -> PipelineResult:
    device = infer_device() if config.device is None else config.device
    LOGGER.info("Pipeline start: device=%s", device)

    out = Path(config.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    result = PipelineResult()

    if config.data_path:
        LOGGER.info("Stage 0: IO — platform=%s", config.platform or "auto-detect")
        adata0 = run_stage0(config)
        result.stage0_output = str(out / "stage0.h5ad")
        adata0.write_h5ad(result.stage0_output)
        LOGGER.info("Stage 0 output: %s", result.stage0_output)
    else:
        result.notes.append("Stage 0 skipped (no --data-path)")

    if not config.skip_alignment and config.data_path:
        LOGGER.info("Stage 1: Alignment")
        try:
            r1 = run_stage1(config, None)
            result.stage1_output = r1.outdir if hasattr(r1, "outdir") else str(out / "aligned.h5ad")
        except Exception as exc:
            result.notes.append(f"Stage 1 failed: {exc}")
            LOGGER.warning("Stage 1 failed: %s", exc)

    if not config.skip_reconstruction:
        LOGGER.info("Stage 2: Reconstruction (model=%s)", config.model_type)
        try:
            result.stage2_output = run_stage2(config, None)
        except Exception as exc:
            result.notes.append(f"Stage 2 failed: {exc}")
            LOGGER.warning("Stage 2 failed: %s", exc)

    if not config.skip_morphogenesis:
        LOGGER.info("Stage 3: Morphogenesis")
        try:
            result.stage3_output = run_stage3(config, None)
        except Exception as exc:
            result.notes.append(f"Stage 3 failed: {exc}")
            LOGGER.warning("Stage 3 failed: %s", exc)

    result_path = out / "pipeline_result.json"
    result_path.write_text(json.dumps(asdict(result), indent=2, default=str))
    LOGGER.info("Pipeline result: %s", result_path)
    return result


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Full Spateo pipeline orchestrator")
    p.add_argument("--data-path", type=str, default=None)
    p.add_argument("--platform", type=str, default=None)
    p.add_argument("--auto-detect", action="store_true")
    p.add_argument("--slices-dir", type=str, default=None)
    p.add_argument("--stage2-path", type=str, default=None)
    p.add_argument("--out-dir", type=str, default="./output/")
    p.add_argument("--skip-alignment", action="store_true")
    p.add_argument("--skip-reconstruction", action="store_true")
    p.add_argument("--skip-morphogenesis", action="store_true")
    p.add_argument("--model-type", type=str, default="surface")
    p.add_argument("--interpolation-method", type=str, default="vtk")
    p.add_argument("--groupby", type=str, default="tissue_type")
    p.add_argument("--device", type=str, default=None)
    p.add_argument("--verbose", action="store_true")
    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_arg_parser().parse_args(argv)
    configure_logging(args.verbose)

    cfg = PipelineConfig(
        data_path=args.data_path,
        platform=args.platform,
        auto_detect=args.auto_detect,
        slices_dir=args.slices_dir,
        stage2_path=args.stage2_path,
        out_dir=args.out_dir,
        skip_alignment=args.skip_alignment,
        skip_reconstruction=args.skip_reconstruction,
        skip_morphogenesis=args.skip_morphogenesis,
        model_type=args.model_type,
        interpolation_method=args.interpolation_method,
        groupby=args.groupby,
        device=args.device,
        verbose=args.verbose,
    )
    result = run_full_pipeline(cfg)
    print(json.dumps(asdict(result), indent=2, default=str))


if __name__ == "__main__":
    main()
