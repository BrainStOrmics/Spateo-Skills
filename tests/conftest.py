"""Shared pytest fixtures for Spateo-Skills tests."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

SKILLS_ROOT = Path(__file__).resolve().parent.parent / "skills"


@pytest.fixture
def skills_dir():
    """Return the skills/ directory path."""
    return SKILLS_ROOT


@pytest.fixture
def syspath_with_skills():
    """Add skills/ to sys.path for import tests."""
    p = str(SKILLS_ROOT)
    if p not in sys.path:
        sys.path.insert(0, p)
    yield
    if p in sys.path:
        sys.path.remove(p)
