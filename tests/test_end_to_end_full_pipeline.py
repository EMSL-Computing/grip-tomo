"""Comprehensive end-to-end smoke test of synthetic workflow.

Pipeline covered:
- Parse `try_these_conditions.txt` via `prepare_multiple_conditions.prepare_many_inputs`.
- Enter one generated condition folder, place a noise MRC file.
- Run `at_single_folder.main` with stubs for merge, phaseflip, reconstruction, filter, and invert steps.

Artifacts validated:
- Condition folder naming pattern from prepare step.
- Merge produced initial stack `.mrcs` file.
- Phase flip produced `_ctfcorrected.mrcs` output.
- Reconstruction produced a `.rec` file retained (ctfcorrected tagged) and moved into reconstruction directory.
- Filter and invert scripts invoked.
- Angle file, defocus file, CSV listing moved to reconstruction directory.
- Original noise MRC archived under `mrc/`.

Side effects are simulated to avoid heavy compute.
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest

# Dynamic paths
UTILS_ROOT = Path(__file__).resolve().parents[1] / "utils"
PREPARE_PATH = UTILS_ROOT / "prepare_multiple_conditions.py"
SINGLE_FOLDER_PATH = UTILS_ROOT / "at_single_folder.py"


def _load(path: Path) -> ModuleType:
    name = path.stem + "_" + os.urandom(4).hex()
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="function")
def prepare_module():
    module = _load(PREPARE_PATH)
    try:
        yield module
    finally:
        sys.modules.pop(module.__name__, None)


@pytest.fixture(scope="function")
def single_folder_module():
    module = _load(SINGLE_FOLDER_PATH)
    try:
        yield module
    finally:
        sys.modules.pop(module.__name__, None)


def test_full_pipeline(tmp_path, monkeypatch, prepare_module, single_folder_module):
    # 1. Create minimal try_these_conditions.txt with one value each
    conds_file = tmp_path / "try_these_conditions.txt"
    conds_file.write_text(
        "doseval 1\nframeval 1\nthickval 10\nphaseval 0\ndefocusval 100\n",
        encoding="utf-8",
    )

    # Template input directory with required files
    template_root = tmp_path / "template_input"
    template_root.mkdir()
    (template_root / "template_simulate-tilt-noise.sh").write_text(
        "#!/bin/bash\necho replace_here\n", encoding="utf-8"
    )
    (template_root / "template_runscript.py").write_text(
        "print('replace_this doseval')\nprint('replace_this frameval')\nprint('replace_this thickval')\nprint('replace_this phaseval')\nprint('replace_this defocusval')\n",
        encoding="utf-8",
    )

    # PDB input directory
    pdb_dir = tmp_path / "pdb_inputs"
    pdb_dir.mkdir()
    (pdb_dir / "proteinA.pdb").write_text("ATOM\n", encoding="utf-8")

    monkeypatch.setenv("GRIPTOMO_TEMPLATE_ROOT", str(template_root))
    monkeypatch.chdir(tmp_path)

    # Run prepare step
    assert prepare_module.prepare_many_inputs(str(conds_file), str(pdb_dir)) is True

    # Determine generated condition folder
    each_input_root = tmp_path / "each_input_folders" / "proteinA"
    # Expected folder name pattern: dose x frame e_thick t_phase p_defocus d
    condition_dirs = [
        d
        for d in each_input_root.iterdir()
        if d.is_dir() and d.name != "template_input"
    ]
    assert len(condition_dirs) == 1
    condition_dir = condition_dirs[0]
    assert condition_dir.name.endswith("d") and "x" in condition_dir.name

    # Place a synthetic noise MRC for angle 0deg
    noise_mrc = condition_dir / "sample_noise_0deg.mrc"
    noise_mrc.write_text("density", encoding="utf-8")

    # 2. Prepare environment overrides for reconstruction pipeline
    merge_script = tmp_path / "merge_mrc_s_in_this_list_noise.py"
    merge_script.write_text("# merge stub\n", encoding="utf-8")
    reconstruct_script = tmp_path / "2_reconstruct_simulated.bash"
    reconstruct_script.write_text("#!/bin/bash\n# reconstruct stub\n", encoding="utf-8")
    phase_template = tmp_path / "in_nm.defocus_template"
    phase_template.write_text(
        "# header\n12 values line\n1 2 3 4 5 6\n", encoding="utf-8"
    )

    monkeypatch.setenv("GRIPTOMO_MERGE_SCRIPT", str(merge_script))
    monkeypatch.setenv("GRIPTOMO_RECONSTRUCT_SCRIPT", str(reconstruct_script))
    monkeypatch.setenv("GRIPTOMO_PHASE_TEMPLATE", str(phase_template))
    monkeypatch.setenv("IMOD_CTFPHASEFLIP_BIN", "ctfphaseflip")
    monkeypatch.setenv("GRIPTOMO_PYTHON_BIN", "python3")

    # Stub util.util.show_time
    util_pkg = ModuleType("util")
    util_submodule = ModuleType("util.util")
    util_submodule.show_time = (
        lambda label, start, end: f"{label} took {end - start:.2f}s"
    )
    util_pkg.util = util_submodule
    monkeypatch.setitem(sys.modules, "util", util_pkg)
    monkeypatch.setitem(sys.modules, "util.util", util_submodule)

    calls: list[tuple[str, ...]] = []

    def fake_run_command(cmd, cwd=None):  # mimic merge/phase/reconstruct side effects
        tcmd = tuple(map(str, cmd))
        calls.append(tcmd)
        cwd_path = Path(cwd) if cwd else condition_dir
        if str(merge_script) in tcmd:
            # Produce initial stack
            (cwd_path / "sample_0p_0deg.mrcs").write_text("stack", encoding="utf-8")
        elif "ctfphaseflip" in tcmd[0]:
            # Produce ctfcorrected variant
            (cwd_path / "sample_0p_0deg_degPhase_0.0_ctfcorrected.mrcs").write_text(
                "ctf", encoding="utf-8"
            )
        elif str(reconstruct_script) in tcmd:
            # Produce reconstruction output
            (cwd_path / "volume_ctfcorrected.rec").write_text("rec", encoding="utf-8")
        return SimpleNamespace(stdout="", stderr="", returncode=0)

    monkeypatch.setattr(single_folder_module, "run_command", fake_run_command)

    filter_calls: list[tuple[Path, float]] = []

    def fake_apply_filters(
        rec_dir, low_pass_level, *, filter_module=None, invert_module=None
    ):
        filter_calls.append((Path(rec_dir), float(low_pass_level)))
        rec_dir = Path(rec_dir)
        rec_dir.mkdir(exist_ok=True)
        (rec_dir / "filtered.mrc").write_text("filtered", encoding="utf-8")
        (rec_dir / "inverted.mrc").write_text("inverted", encoding="utf-8")

    monkeypatch.setattr(single_folder_module, "apply_filters", fake_apply_filters)

    # 3. Run single folder pipeline
    monkeypatch.chdir(condition_dir)
    monkeypatch.setattr(
        sys, "argv", [str(SINGLE_FOLDER_PATH), "0", "0", "1", "0.5"], raising=False
    )
    single_folder_module.main()

    # reconstruction directory name logic
    rec_dir = condition_dir / "0_0_1.0"
    assert rec_dir.is_dir(), "Reconstruction directory missing"

    # Validate moved outputs
    assert (rec_dir / "df_of_all_mrc_files_sorted_by_angle.csv").exists()
    assert (rec_dir / "angles.rawtlt").exists()
    assert (rec_dir / "in_nm.defocus").exists()
    assert (rec_dir / "sample_0p_0deg.mrcs").exists(), "Merged stack missing"
    assert (rec_dir / "sample_0p_0deg_degPhase_0.0_ctfcorrected.mrcs").exists(), (
        "CTF-corrected stack missing"
    )
    assert (rec_dir / "volume_ctfcorrected_0_0_1.0.rec").exists(), (
        "Renamed .rec file missing"
    )
    assert (rec_dir / "filtered.mrc").exists(), "Filtered output missing"
    assert (rec_dir / "inverted.mrc").exists(), "Inverted output missing"

    # Archived original mrc
    archive_mrc_dir = condition_dir / "mrc"
    assert archive_mrc_dir.is_dir() and list(archive_mrc_dir.glob("*.mrc")), (
        "Noise MRC not archived"
    )

    # Commands invoked
    flat = [" ".join(c) for c in calls]
    assert any("merge_mrc_s_in_this_list_noise.py" in c for c in flat)
    assert any(c.startswith("ctfphaseflip") for c in flat)
    assert any("2_reconstruct_simulated.bash" in c for c in flat)
    assert filter_calls and filter_calls[0][0] == rec_dir
