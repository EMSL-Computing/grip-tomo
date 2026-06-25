"""Tests for run_at_many_input_folders_that_has_both_many_pdb_and_many_conditions."""

from __future__ import annotations

import importlib.util
import os
import sys
import uuid
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    REPO_ROOT
    / "utils"
    / "run_at_many_input_folders_that_has_both_many_pdb_and_many_conditions.py"
)


def _load_module(monkeypatch, workdir: Path):
    """Import the legacy run_at_many_input_folders module for testing."""
    module_name = f"run_at_many_input_folders_{uuid.uuid4().hex}"
    spec = importlib.util.spec_from_file_location(module_name, MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load run_at_many_input_folders module spec")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    module.starting_dir = str(workdir)
    return module


def test_run_at_many_folders_organizes_outputs(tmp_path, monkeypatch):
    """Ensure run_at_many_folders creates output/input directories and moves files."""
    each_input = tmp_path / "each_input_folders"
    condition_dir = each_input / "pdb_one" / "condition_a"
    condition_dir.mkdir(parents=True)

    (condition_dir / "cisTEM.log").write_text("log\n")
    (condition_dir / "simulate_1x10e_3t_4p").write_text("sim\n")
    (condition_dir / "misc.txt").write_text("misc\n")

    module = _load_module(monkeypatch, each_input)

    original_cwd = os.getcwd()
    try:
        module.run_at_many_folders("unused.py")
    finally:
        os.chdir(original_cwd)

    assert sorted(os.listdir(condition_dir)) == ["input", "output"]
    output_files = sorted(os.listdir(condition_dir / "output"))
    assert output_files == ["cisTEM.log", "simulate_1x10e_3t_4p"]
    input_files = sorted(os.listdir(condition_dir / "input"))
    assert input_files == ["misc.txt"]


def test_run_at_many_folders_skips_existing_output(tmp_path, monkeypatch):
    """Ensure run_at_many_folders leaves folders untouched when output exists."""
    each_input = tmp_path / "each_input_folders"
    condition_dir = each_input / "pdb_one" / "condition_a"
    (condition_dir / "output").mkdir(parents=True)

    (condition_dir / "output" / "existing.txt").write_text("old\n")
    (condition_dir / "cisTEM.log").write_text("log\n")

    module = _load_module(monkeypatch, each_input)

    original_cwd = os.getcwd()
    try:
        module.run_at_many_folders("unused.py")
    finally:
        os.chdir(original_cwd)

    assert (condition_dir / "output" / "existing.txt").exists()
    assert not (condition_dir / "input").exists()
    assert (condition_dir / "cisTEM.log").exists()
