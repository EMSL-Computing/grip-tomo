"""Tests for utils.use_parallelized_cores helper."""

from __future__ import annotations

import importlib.util
import os
import sys
import uuid
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "utils" / "use_parallelized_cores.py"


def _load_module(monkeypatch, workdir: Path):
    """Import use_parallelized_cores into an isolated module namespace."""
    module_name = f"use_parallelized_cores_{uuid.uuid4().hex}"
    spec = importlib.util.spec_from_file_location(module_name, MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load use_parallelized_cores module spec")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    module.starting_dir = str(workdir)
    return module


def test_modified_run_submits_all_simulations(tmp_path, monkeypatch):
    """Ensure modified_run_at_many_folders submits each simulate subdirectory."""
    workdir = tmp_path / "each_input_folders"
    output_dir = workdir / "pdb_one" / "condition_a" / "output"
    (output_dir / "simulate_a").mkdir(parents=True)
    (output_dir / "simulate_b").mkdir()

    module = _load_module(monkeypatch, workdir)

    submitted_calls: list[tuple] = []

    def fake_process(
        directory, run_this_py, min_tilt, max_tilt, increment, low_pass_level
    ):
        submitted_calls.append(
            (directory, run_this_py, min_tilt, max_tilt, increment, low_pass_level)
        )
        return None

    monkeypatch.setattr(module, "process_directory", fake_process)

    class _ImmediateFuture:
        def __init__(self, fn, args, kwargs):
            self._result = fn(*args, **kwargs)

        def result(self):
            return self._result

    class _FakeExecutor:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def submit(self, fn, *args, **kwargs):
            return _ImmediateFuture(fn, args, kwargs)

    monkeypatch.setattr(
        module.concurrent.futures,
        "ProcessPoolExecutor",
        lambda *a, **k: _FakeExecutor(),
    )
    monkeypatch.setattr(
        module.concurrent.futures, "as_completed", lambda futures: futures
    )

    original_cwd = os.getcwd()
    try:
        module.modified_run_at_many_folders("run_at_single_folder.py", -60, 60, 3, 0.16)
    finally:
        os.chdir(original_cwd)

    directories = {Path(call[0]).name for call in submitted_calls}
    assert directories == {"simulate_a", "simulate_b"}


def test_modified_run_exits_without_output(tmp_path, monkeypatch):
    """Ensure modified_run_at_many_folders exits when output is missing."""
    workdir = tmp_path / "each_input_folders"
    (workdir / "pdb_one" / "condition_a").mkdir(parents=True)

    module = _load_module(monkeypatch, workdir)

    class _FakeExecutor:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def submit(self, fn, *args, **kwargs):
            raise AssertionError("submit should not be called when output is missing")

    monkeypatch.setattr(
        module.concurrent.futures,
        "ProcessPoolExecutor",
        lambda *a, **k: _FakeExecutor(),
    )

    with pytest.raises(SystemExit) as exc_info:
        module.modified_run_at_many_folders("run.py", -60, 60, 3, 0.16)

    assert exc_info.value.code == 1


def test_process_directory_invokes_python(tmp_path, monkeypatch):
    """Ensure process_directory shells out to the expected helper script."""
    module = _load_module(monkeypatch, tmp_path)

    simulate_dir = tmp_path / "simulate_dir"
    simulate_dir.mkdir()

    captured: dict[str, object] = {}

    def fake_run(cmd, check, cwd):
        captured["cmd"] = cmd
        captured["check"] = check
        captured["cwd"] = cwd

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    original_cwd = os.getcwd()
    try:
        module.process_directory(
            str(simulate_dir),
            "run_at_single_folder.py",
            -60,
            60,
            3,
            0.16,
        )
    finally:
        os.chdir(original_cwd)

    cmd_list = captured["cmd"]
    assert isinstance(cmd_list, list)
    assert cmd_list[0] == sys.executable
    assert cmd_list[1] == "run_at_single_folder.py"
    assert cmd_list[2:] == ["-60", "60", "3", "0.16"]
    assert captured["check"] is False
    assert captured["cwd"] == str(simulate_dir)


def test_process_directory_respects_python_env(tmp_path, monkeypatch):
    """Ensure process_directory honors GRIPTOMO_PYTHON_BIN override."""
    module = _load_module(monkeypatch, tmp_path)

    simulate_dir = tmp_path / "simulate_dir"
    simulate_dir.mkdir()

    monkeypatch.setenv("GRIPTOMO_PYTHON_BIN", "/custom/python")

    captured: dict[str, object] = {}

    def fake_run(cmd, check, cwd):
        captured["cmd"] = cmd

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    module.process_directory(
        str(simulate_dir),
        "run_at_single_folder.py",
        -60,
        60,
        3,
        0.16,
    )

    assert captured["cmd"][0] == "/custom/python"


def test_modified_run_honors_max_workers_env(tmp_path, monkeypatch):
    """Ensure modified_run_at_many_folders forwards GRIPTOMO_MAX_WORKERS."""
    workdir = tmp_path / "each_input_folders"
    output_dir = workdir / "pdb_one" / "condition_a" / "output"
    (output_dir / "simulate_a").mkdir(parents=True)

    module = _load_module(monkeypatch, workdir)
    monkeypatch.setenv("GRIPTOMO_MAX_WORKERS", "4")

    class _ImmediateFuture:
        def __init__(self, fn, args, kwargs):
            self._result = None

        def result(self):
            return self._result

    captured_kwargs: dict[str, object] = {}

    class _FakeExecutor:
        def __init__(self, **kwargs):
            captured_kwargs.update(kwargs)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def submit(self, fn, *args, **kwargs):
            return _ImmediateFuture(fn, args, kwargs)

    monkeypatch.setattr(
        module.concurrent.futures,
        "ProcessPoolExecutor",
        lambda **kwargs: _FakeExecutor(**kwargs),
    )
    monkeypatch.setattr(
        module.concurrent.futures, "as_completed", lambda futures: futures
    )
    monkeypatch.setattr(module, "process_directory", lambda *a, **k: None)

    module.modified_run_at_many_folders("run.py", -60, 60, 3, 0.16)

    assert captured_kwargs.get("max_workers") == 4
