"""Lazy import helpers for Spateo and heavy scientific dependencies."""

from __future__ import annotations

from typing import Any


def _import_spateo() -> Any:
    """Import spateo lazily so syntax checks don't require installation."""
    try:
        import spateo as st  # type: ignore
        return st
    except ImportError as exc:
        raise ImportError(
            "spateo is required. Install from envs/requirements.txt or "
            "Spateo-Skills/00_env_setup/."
        ) from exc


def _import_spateo_plus() -> Any | None:
    """Import spateo_plus if available; return None otherwise."""
    try:
        import spateo_plus as sp  # type: ignore
        return sp
    except Exception:
        return None


def _import_pyvista() -> Any:
    """Import pyvista lazily."""
    try:
        import pyvista as pv  # type: ignore
        return pv
    except ImportError as exc:
        raise ImportError(
            "pyvista is required for 3D model operations."
        ) from exc
