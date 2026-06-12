"""Verify all skill modules are syntactically valid Python.

These tests check that modules parse and compile without Spateo installed.
They do NOT test runtime behavior (that requires the full Spateo env).
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

SKILLS_ROOT = Path(__file__).resolve().parent.parent / "skills"


def _skill_modules():
    """Yield all .py files under skills/ (excluding __init__.py)."""
    for py in SKILLS_ROOT.rglob("*.py"):
        if py.name == "__init__.py":
            continue
        yield py


@pytest.mark.parametrize("module", list(_skill_modules()), ids=lambda p: p.relative_to(SKILLS_ROOT))
def test_syntax(module: Path):
    """Every skill module must be valid Python syntax."""
    source = module.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(module))
    # Verify at least one class or function definition exists
    has_definition = any(
        isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
        for node in ast.walk(tree)
    )
    assert has_definition, f"{module.relative_to(SKILLS_ROOT)} has no class/function definitions"


def test_shared_imports_compile():
    """00_shared modules must compile (no Spateo needed for lazy import helpers)."""
    shared = SKILLS_ROOT / "00_shared"
    for mod in ["lazy_imports.py", "data_helpers.py"]:
        source = (shared / mod).read_text(encoding="utf-8")
        ast.parse(source)


def test_env_setup_files_exist():
    """Environment setup files must exist."""
    env = SKILLS_ROOT / "00_env_setup"
    assert (env / "setup.sh").exists(), "setup.sh missing"
    assert (env / "environment.yml").exists(), "environment.yml missing"
    assert (env / "uv_requirements.txt").exists(), "uv_requirements.txt missing"
    assert (env / "uv_pyproject.toml").exists(), "uv_pyproject.toml missing"
    assert (env / "skill.md").exists(), "skill.md missing"


def test_skill_md_files_exist():
    """Each skill directory should have a SKILL.md or skill.md."""
    for skill_dir in SKILLS_ROOT.iterdir():
        if skill_dir.is_dir() and not skill_dir.name.startswith("."):
            has_any = (skill_dir / "SKILL.md").exists() or (skill_dir / "skill.md").exists()
            assert has_any, f"{skill_dir.name}/ missing SKILL.md or skill.md"
