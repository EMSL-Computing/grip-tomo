"""Tests for single-condition run_at_single_input_folder helper."""

from __future__ import annotations

import importlib.util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "utils" / "run_at_single_input_folder.py"


def _load_module(module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)  # type: ignore[call-arg]
    return module


def test_run_at_single_moves_outputs(tmp_path, monkeypatch):
    monkeypatch.setenv("GRIPTOMO_OUTPUT_DIR", str(tmp_path / "out"))
    """Ensure run_at_single_input_folder.move_outputs relocates expected files."""
    module = _load_module("run_at_single_input_folder_module")

    workdir = tmp_path / "input"
    workdir.mkdir()
    (workdir / "cistem.log").write_text("log", encoding="utf-8")
    (workdir / "simulate_1x1e_1t_1p").write_text("sim", encoding="utf-8")
    (workdir / "keep.txt").write_text("keep", encoding="utf-8")

    module.move_outputs(workdir, module.OUTPUT_DIR)

    output_dir = Path(module.OUTPUT_DIR)
    assert sorted(p.name for p in output_dir.iterdir()) == [
        "cistem.log",
        "simulate_1x1e_1t_1p",
    ]
    assert (workdir / "keep.txt").exists()


def test_move_outputs_uses_default_output(tmp_path, monkeypatch):
    monkeypatch.delenv("GRIPTOMO_OUTPUT_DIR", raising=False)
    """Verify default OUTPUT_DIR (../output) is used when env is unset."""
    module = _load_module("run_at_single_input_folder_default")

    workdir = tmp_path / "input"
    workdir.mkdir()
    (workdir / "cistem.err").write_text("err", encoding="utf-8")

    module.move_outputs(workdir, module.OUTPUT_DIR)
    expected_dir = (workdir / module.OUTPUT_DIR).resolve()
    assert module.OUTPUT_DIR == Path("../output")
    assert (expected_dir / "cistem.err").exists()
