"""Smoke tests for utils.at_single_folder."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

MODULE_NAME = "utils.at_single_folder"


def _load_module():
    module = importlib.import_module(MODULE_NAME)
    return importlib.reload(module)


@pytest.fixture
def at_single_folder_module():
    module = _load_module()
    try:
        yield module
    finally:
        sys.modules.pop(module.__name__, None)


def test_generate_angles_inclusive(at_single_folder_module):
    module = at_single_folder_module
    angles = module.generate_angles(-1, 1, 1)
    assert angles == [-1.0, 0.0, 1.0]


def test_build_pipeline_paths_env_override(
    at_single_folder_module, tmp_path, monkeypatch
):
    module = at_single_folder_module
    defaults = module.build_pipeline_paths(Path(module.__file__).resolve())

    merge_override = tmp_path / "custom_merge.py"
    merge_override.write_text("", encoding="utf-8")
    monkeypatch.setenv("GRIPTOMO_MERGE_SCRIPT", str(merge_override))

    paths = module.build_pipeline_paths(Path(module.__file__).resolve())

    assert paths.merge_script == merge_override
    assert paths.reconstruct_script == defaults.reconstruct_script


def test_main_smoke_execution(at_single_folder_module, tmp_path, monkeypatch, capsys):
    module = at_single_folder_module

    util_pkg = ModuleType("util")
    util_submodule = ModuleType("util.util")

    def show_time(label, start, end):
        return f"{label} took {end - start:.2f}s"

    util_submodule.show_time = show_time
    util_pkg.util = util_submodule

    monkeypatch.setitem(sys.modules, "util", util_pkg)
    monkeypatch.setitem(sys.modules, "util.util", util_submodule)

    merge_script = tmp_path / "merge.py"
    phase_template = tmp_path / "phase_template.txt"
    reconstruct_script = tmp_path / "reconstruct.sh"
    merge_script.write_text("", encoding="utf-8")
    phase_template.write_text(
        "# template\n12 header line\n1 2 3 4 5 6\n",
        encoding="utf-8",
    )
    reconstruct_script.write_text("", encoding="utf-8")

    monkeypatch.setenv("GRIPTOMO_MERGE_SCRIPT", str(merge_script))
    monkeypatch.setenv("GRIPTOMO_PHASE_TEMPLATE", str(phase_template))
    monkeypatch.setenv("GRIPTOMO_RECONSTRUCT_SCRIPT", str(reconstruct_script))
    monkeypatch.setenv("IMOD_CTFPHASEFLIP_BIN", "ctfphaseflip")
    monkeypatch.setenv("GRIPTOMO_PYTHON_BIN", "python3")

    import utils.filter_mrc as filter_mrc_mod
    import utils.invert_mrc_density as invert_mrc_mod

    def fake_filter(level, *, command_runner=None):
        Path("filtered.mrc").write_text("filtered", encoding="utf-8")

    def fake_invert(directory, *, command_runner=None):
        Path(directory, "inverted.mrc").write_text("inverted", encoding="utf-8")

    monkeypatch.setattr(filter_mrc_mod, "filter_mrc_files", fake_filter)
    monkeypatch.setattr(invert_mrc_mod, "invert_directory", fake_invert)

    commands: list[tuple[tuple[str, ...], Path]] = []

    def fake_run_command(cmd, cwd=None):
        path_cwd = Path(cwd) if cwd else tmp_path
        commands.append((tuple(cmd), path_cwd))
        if str(merge_script) in cmd:
            (path_cwd / "sample_0p_0deg.mrcs").write_text("", encoding="utf-8")
        return SimpleNamespace(stdout="", stderr="", returncode=0)

    monkeypatch.setattr(module, "run_command", fake_run_command)

    monkeypatch.chdir(tmp_path)
    (tmp_path / "sample_noise_0deg.mrc").write_text("data", encoding="utf-8")

    monkeypatch.setattr(
        sys, "argv", [str(Path(module.__file__)), "0", "0", "1", "0.5"], raising=False
    )

    module.main()

    rec_dir_name = "0_0_1.0"
    rec_dir = tmp_path / rec_dir_name
    assert rec_dir.is_dir()
    assert (rec_dir / "df_of_all_mrc_files_sorted_by_angle.csv").exists()
    assert (rec_dir / "angles.rawtlt").exists()
    assert (rec_dir / "in_nm.defocus").exists()
    assert (rec_dir / "sample_0p_0deg.mrcs").exists()

    archive_dir = tmp_path / "mrc"
    assert (archive_dir / "sample_noise_0deg.mrc").exists()

    assert any(cmd[0][0] == "ctfphaseflip" for cmd in commands)
    assert any(str(merge_script) in cmd[0] for cmd in commands)

    out = capsys.readouterr().out
    assert "reconstruction took" in out


def test_phase_template_env_missing_uses_default(
    at_single_folder_module, tmp_path, monkeypatch
):
    """If GRIPTOMO_PHASE_TEMPLATE is unset, build_pipeline_paths should fall back to default template path."""
    module = at_single_folder_module
    # Ensure variable absent
    monkeypatch.delenv("GRIPTOMO_PHASE_TEMPLATE", raising=False)
    defaults = module.build_pipeline_paths(Path(module.__file__).resolve())
    # Provide other required overrides to avoid unrelated failures
    merge_override = tmp_path / "merge.py"
    reconstruct_override = tmp_path / "reconstruct.sh"
    for p in [merge_override, reconstruct_override]:
        p.write_text("", encoding="utf-8")
    monkeypatch.setenv("GRIPTOMO_MERGE_SCRIPT", str(merge_override))
    monkeypatch.setenv("GRIPTOMO_RECONSTRUCT_SCRIPT", str(reconstruct_override))
    # Validate fallback
    paths = module.build_pipeline_paths(Path(module.__file__).resolve())
    assert paths.phase_template == defaults.phase_template
