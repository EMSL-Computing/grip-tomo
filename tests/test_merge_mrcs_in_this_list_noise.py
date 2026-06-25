"""Tests for merge_mrc_s_in_this_list_noise helper functions."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import mrcfile
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "utils" / "merge_mrc_s_in_this_list_noise.py"


def _load_module(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)  # type: ignore[call-arg]
    return module


def test_check_filters_zero_mean(tmp_path, monkeypatch):
    module = _load_module("merge_mrc_module")
    monkeypatch.chdir(tmp_path)

    good_path = tmp_path / "tilt_good.mrc"
    bad_path = tmp_path / "tilt_bad.mrc"

    with mrcfile.new(good_path, overwrite=True) as mrc:
        mrc.set_data(np.ones((1, 2, 2), dtype=np.float32))

    with mrcfile.new(bad_path, overwrite=True) as mrc:
        mrc.set_data(np.zeros((1, 2, 2), dtype=np.float32))

    df = pd.DataFrame(
        {
            "mrc_filename": [str(good_path), str(bad_path)],
            "angle": [0.0, 10.0],
        }
    )

    result = module.check_which_mrc_file_can_be_merged_by_averaging(df)
    assert list(result["mrc_filename"]) == [str(good_path)]
    csv_path = Path("df_of_proper_mrc_files_sorted_by_angle.csv")
    assert csv_path.exists()


def test_when_all_mrc_files_builds_command(tmp_path, monkeypatch):
    module = _load_module("merge_mrc_module_cmd")
    monkeypatch.chdir(tmp_path)

    mrc_path = tmp_path / "tilt_single.mrc"
    with mrcfile.new(mrc_path, overwrite=True) as mrc:
        mrc.set_data(np.ones((1, 2, 2), dtype=np.float32))

    df = pd.DataFrame({"mrc_filename": [str(mrc_path)], "angle": [0.0]})

    captured = {}

    def fake_system(cmd):
        captured["command"] = cmd
        return 0

    monkeypatch.setattr(module.os, "system", fake_system)

    module.when_all_mrc_files_are_proper_to_be_merged(df, "/opt/eman2/bin/e2proc2d.py")
    assert "e2proc2d.py" in captured["command"]
    assert str(mrc_path) in captured["command"]
    assert captured["command"].strip().endswith(".mrcs")


def test_locate_e2proc2d_respects_env(tmp_path, monkeypatch):
    module = _load_module("merge_mrc_module_locate")
    fake_bin = tmp_path / "e2proc2d.py"
    fake_bin.write_text("#!/bin/sh\n", encoding="utf-8")
    fake_bin.chmod(0o755)

    monkeypatch.setenv("EMAN2_E2PROC2D", str(fake_bin))
    assert module._locate_e2proc2d() == str(fake_bin)
