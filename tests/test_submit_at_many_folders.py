"""Tests for utils.submit_at_many_folders helper."""

from __future__ import annotations

import importlib.util
import sys
import types
import uuid
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "utils" / "submit_at_many_folders.py"


def _load_submit_module(monkeypatch, workdir: Path):
    """Import submit_at_many_folders with a controlled working directory."""
    dummy_pretty_errors = types.ModuleType("pretty_errors")
    monkeypatch.setitem(sys.modules, "pretty_errors", dummy_pretty_errors)
    monkeypatch.chdir(workdir)

    module_name = f"submit_at_many_folders_{uuid.uuid4().hex}"
    spec = importlib.util.spec_from_file_location(module_name, MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load submit_at_many_folders module spec")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def test_submit_at_many_folders_invokes_sbatch(tmp_path, monkeypatch):
    """Ensure submit_at_many_folders issues sbatch for each condition directory."""
    workdir = tmp_path / "each_input_folders"
    workdir.mkdir()

    (workdir / "pdb_one" / "condition_a").mkdir(parents=True)
    (workdir / "pdb_one" / "template_input").mkdir()
    (workdir / "pdb_two" / "condition_b").mkdir(parents=True)

    module = _load_submit_module(monkeypatch, workdir)

    calls: list[tuple[str, str]] = []

    original_system = module.os.system

    def fake_system(cmd: str) -> int:
        calls.append((module.os.getcwd(), cmd))
        return 0

    module.os.system = fake_system  # type: ignore[assignment]
    start_cwd = module.os.getcwd()

    try:
        module.submit_at_many_folders()
        assert module.os.getcwd() == start_cwd
    finally:
        module.os.system = original_system  # type: ignore[assignment]
        module.os.chdir(start_cwd)

    assert {command for _, command in calls} == {"sbatch run.sbatch"}
    executed_dirs = {cwd for cwd, _ in calls}
    expected_dirs = {
        str(workdir / "pdb_one" / "condition_a"),
        str(workdir / "pdb_two" / "condition_b"),
    }
    assert executed_dirs == expected_dirs


def test_submit_at_many_folders_respects_env_override(tmp_path, monkeypatch):
    """Ensure submit_at_many_folders honors SLURM_SBATCH override."""
    workdir = tmp_path / "each_input_folders"
    (workdir / "pdb_one" / "condition_a").mkdir(parents=True)
    module = _load_submit_module(monkeypatch, workdir)

    monkeypatch.setenv("SLURM_SBATCH", "custom_sbatch")

    commands: list[str] = []
    original_system = module.os.system
    start_cwd = module.os.getcwd()

    def fake_system(cmd: str) -> int:
        commands.append(cmd)
        return 0

    module.os.system = fake_system  # type: ignore[assignment]

    try:
        module.submit_at_many_folders()
        assert module.os.getcwd() == start_cwd
    finally:
        module.os.system = original_system  # type: ignore[assignment]
        module.os.chdir(start_cwd)

    assert commands == ["custom_sbatch run.sbatch"]


def test_submit_at_many_folders_stops_on_stars(tmp_path, monkeypatch):
    """Ensure submit_at_many_folders exits when stars directory is discovered."""
    workdir = tmp_path / "each_input_folders"
    (workdir / "pdb_one" / "stars").mkdir(parents=True)
    module = _load_submit_module(monkeypatch, workdir)

    start_cwd = module.os.getcwd()
    try:
        with pytest.raises(SystemExit) as exc_info:
            module.submit_at_many_folders()
    finally:
        module.os.chdir(start_cwd)

    assert exc_info.value.code == 1
