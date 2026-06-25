"""Tests for utils.prepare_multiple_conditions helper functions."""

from __future__ import annotations

import importlib.util
import sys
import uuid
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "utils" / "prepare_multiple_conditions.py"


def _load_prepare_multiple_conditions():
    """Load the legacy prepare_multiple_conditions module for testing."""
    module_name = f"prepare_multiple_conditions_{uuid.uuid4().hex}"
    spec = importlib.util.spec_from_file_location(module_name, MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load prepare_multiple_conditions module spec")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def test_prepare_many_inputs_creates_condition_directories(tmp_path, monkeypatch):
    """Ensure prepare_many_inputs copies template and materializes condition grids."""
    module = _load_prepare_multiple_conditions()

    template_root = tmp_path / "template_source"
    template_input = template_root / "template_input"
    template_input.mkdir(parents=True)
    (template_input / "template_simulate-tilt-noise.sh").write_text(
        "run replace_here\n"
    )
    (template_input / "template_runscript.py").write_text(
        "set doseval replace_this\n"
        "set frameval replace_this\n"
        "set thickval replace_this\n"
        "set phaseval replace_this\n"
        "set defocusval replace_this\n"
    )

    module.template_input_path = str(template_input)

    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("FORCE", "1")

    stale_dir = tmp_path / "each_input_folders"
    stale_dir.mkdir()
    (stale_dir / "stale.txt").write_text("old\n")

    conditions_file = tmp_path / "conditions.txt"
    conditions_file.write_text(
        "doseval 1 2\nframeval 10\nthickval 3\nphaseval 4\ndefocusval 5\n"
    )

    pdb_dir = tmp_path / "pdb_inputs"
    pdb_dir.mkdir()
    (pdb_dir / "sample.pdb").write_text("ATOM\n")

    result = module.prepare_many_inputs(str(conditions_file), str(pdb_dir))

    assert result is True
    assert not (stale_dir / "stale.txt").exists()

    sample_root = tmp_path / "each_input_folders" / "sample"
    template_copy = sample_root / "template_input"
    assert template_copy.exists()
    simulate_script = template_copy / "simulate-tilt-noise.sh"
    assert simulate_script.exists()
    simulate_content = simulate_script.read_text()
    assert "replace_here" not in simulate_content
    assert "../sample.pdb" in simulate_content

    condition_dirs = {
        path.name: path
        for path in sample_root.iterdir()
        if path.is_dir() and path.name != "template_input"
    }
    assert set(condition_dirs) == {"1x10e_3t_4p_5d", "2x10e_3t_4p_5d"}

    for name, directory in condition_dirs.items():
        runscript_path = directory / "runscript.py"
        assert runscript_path.exists()
        runscript_content = runscript_path.read_text()
        assert "replace_this" not in runscript_content
        assert "set frameval 10" in runscript_content
        assert "set thickval 3" in runscript_content
        assert "set phaseval 4" in runscript_content
        assert "set defocusval 5" in runscript_content
        if name.startswith("1x"):
            assert "set doseval 1" in runscript_content
        elif name.startswith("2x"):
            assert "set doseval 2" in runscript_content
        else:
            pytest.fail(f"Unexpected condition directory: {name}")


def test_prepare_many_inputs_returns_false_without_pdbs(tmp_path, monkeypatch):
    """Ensure prepare_many_inputs gracefully handles missing PDB inputs."""
    module = _load_prepare_multiple_conditions()

    template_input = tmp_path / "template_input"
    template_input.mkdir()
    (template_input / "template_simulate-tilt-noise.sh").write_text(
        "run replace_here\n"
    )
    (template_input / "template_runscript.py").write_text("set doseval replace_this\n")

    module.template_input_path = str(template_input)

    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("FORCE", raising=False)

    conditions_file = tmp_path / "conditions.txt"
    conditions_file.write_text(
        "doseval 1\nframeval 10\nthickval 3\nphaseval 4\ndefocusval 5\n"
    )

    inputs_dir = tmp_path / "inputs"
    inputs_dir.mkdir()

    result = module.prepare_many_inputs(str(conditions_file), str(inputs_dir))

    assert result is False
    output_root = tmp_path / "each_input_folders"
    assert output_root.exists()
    assert not any(output_root.iterdir())


def test_prepare_many_inputs_refuses_overwrite_without_force(tmp_path, monkeypatch):
    """Ensure prepare_many_inputs requires FORCE=1 before clobbering output folders."""
    module = _load_prepare_multiple_conditions()

    template_input = tmp_path / "template_input"
    template_input.mkdir()
    (template_input / "template_simulate-tilt-noise.sh").write_text(
        "run replace_here\n"
    )
    (template_input / "template_runscript.py").write_text("set doseval replace_this\n")

    module.template_input_path = str(template_input)

    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("FORCE", raising=False)

    output_root = tmp_path / "each_input_folders"
    output_root.mkdir()
    (output_root / "stale.txt").write_text("old\n")

    conditions_file = tmp_path / "conditions.txt"
    conditions_file.write_text(
        "doseval 1\nframeval 10\nthickval 3\nphaseval 4\ndefocusval 5\n"
    )

    pdb_dir = tmp_path / "inputs"
    pdb_dir.mkdir()
    (pdb_dir / "sample.pdb").write_text("ATOM\n")

    with pytest.raises(FileExistsError):
        module.prepare_many_inputs(str(conditions_file), str(pdb_dir))

    assert (output_root / "stale.txt").exists()
