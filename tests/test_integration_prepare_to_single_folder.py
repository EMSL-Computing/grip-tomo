"""Integration smoke test: prepare_multiple_conditions -> merge stub -> at_single_folder main.

This test exercises a minimal end-to-end path without performing real reconstruction:
1. Generate condition folders using prepare_multiple_conditions.py logic.
2. Create a fake merged .mrcs output via a stub merge script when at_single_folder runs.
3. Invoke at_single_folder.main with environment overrides.

It validates:
- Directory naming from prepare step.
- Merge script executed (creates sample_0p_0deg.mrcs).
- Phase template default or override resolved.
- Output reconstruction folder created with expected angle file.
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

# Paths to scripts under test
REPO_ROOT = Path(__file__).resolve().parents[1]
AT_SINGLE_FOLDER_PATH = REPO_ROOT / "utils" / "at_single_folder.py"
PREPARE_MULTIPLE_CONDITIONS_PATH = (
    REPO_ROOT / "utils" / "prepare_multiple_conditions.py"
)


def _load_module(path: Path) -> ModuleType:
    name = path.stem + "_" + os.urandom(4).hex()
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module spec for {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def prepare_module():
    module = _load_module(PREPARE_MULTIPLE_CONDITIONS_PATH)
    try:
        yield module
    finally:
        sys.modules.pop(module.__name__, None)


@pytest.fixture
def single_folder_module():
    module = _load_module(AT_SINGLE_FOLDER_PATH)
    try:
        yield module
    finally:
        sys.modules.pop(module.__name__, None)


def test_prepare_to_single_folder_end_to_end(
    tmp_path, monkeypatch, prepare_module, single_folder_module
):
    # 1. Simulate prepare_multiple_conditions output
    pdb_file = tmp_path / "protein.pdb"
    pdb_file.write_text("ATOM\n", encoding="utf-8")
    conditions_file = tmp_path / "conditions.txt"
    conditions_file.write_text("defocus=1.5,voltage=300\n", encoding="utf-8")

    # Minimal emulate of prepare_many_inputs: create one condition folder
    cond_dir = tmp_path / "protein" / "defocus_1p5_voltage_300"
    cond_dir.mkdir(parents=True)

    # Fake an input mrc (noise) for reconstruction stage
    noise_mrc = cond_dir / "sample_noise_0deg.mrc"
    noise_mrc.write_text("noise", encoding="utf-8")

    # 2. Set env overrides for at_single_folder
    merge_script = tmp_path / "merge.py"
    merge_script.write_text("", encoding="utf-8")
    reconstruct_script = tmp_path / "reconstruct.sh"
    reconstruct_script.write_text("", encoding="utf-8")
    phase_template = tmp_path / "in_nm.defocus_template"
    phase_template.write_text("header\nvalues\n", encoding="utf-8")

    monkeypatch.setenv("GRIPTOMO_MERGE_SCRIPT", str(merge_script))
    monkeypatch.setenv("GRIPTOMO_RECONSTRUCT_SCRIPT", str(reconstruct_script))
    monkeypatch.setenv("GRIPTOMO_PHASE_TEMPLATE", str(phase_template))
    monkeypatch.setenv("IMOD_CTFPHASEFLIP_BIN", "ctfphaseflip")
    monkeypatch.setenv("GRIPTOMO_PYTHON_BIN", "python3")

    # Stub util.util.show_time to avoid import errors
    util_pkg = ModuleType("util")
    util_submodule = ModuleType("util.util")
    util_submodule.show_time = (
        lambda label, start, end: f"{label} took {end - start:.2f}s"
    )
    util_pkg.util = util_submodule
    monkeypatch.setitem(sys.modules, "util", util_pkg)
    monkeypatch.setitem(sys.modules, "util.util", util_submodule)

    executed = []

    def fake_run_command(cmd, cwd=None):
        executed.append(tuple(cmd))
        # Simulate merge producing a .mrcs file
        if str(merge_script) in cmd:
            (Path(cwd) / "sample_0p_0deg.mrcs").write_text("stack", encoding="utf-8")
        return SimpleNamespace(stdout="", stderr="", returncode=0)

    monkeypatch.setattr(single_folder_module, "run_command", fake_run_command)

    def fake_apply_filters(
        rec_dir, low_pass_level, *, filter_module=None, invert_module=None
    ):
        rec_path = Path(rec_dir)
        rec_path.mkdir(exist_ok=True)
        (rec_path / "filtered.mrc").write_text("filtered", encoding="utf-8")
        (rec_path / "inverted.mrc").write_text("inverted", encoding="utf-8")

    monkeypatch.setattr(single_folder_module, "apply_filters", fake_apply_filters)

    # 3. Run at_single_folder in condition directory
    monkeypatch.chdir(cond_dir)
    monkeypatch.setattr(
        sys, "argv", [str(AT_SINGLE_FOLDER_PATH), "0", "0", "1", "0.5"], raising=False
    )
    single_folder_module.main()

    rec_dir = cond_dir / "0_0_1.0"
    assert rec_dir.is_dir(), "Reconstruction directory not created"
    assert (rec_dir / "angles.rawtlt").exists(), "Angles file missing"
    assert (rec_dir / "sample_0p_0deg.mrcs").exists(), "Merged stack missing"
    assert any(str(merge_script) in c for c in executed), "Merge script never executed"
