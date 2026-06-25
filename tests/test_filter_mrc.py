import subprocess

import importlib

import mrcfile
import numpy as np


def _load_filter_module():
    return importlib.import_module("utils.filter_mrc")


def test_filter_mrc_uses_imod_clip_override(tmp_path, monkeypatch):
    """Ensure IMOD clip env override is honored when filtering MRCs."""
    module = _load_filter_module()
    monkeypatch.chdir(tmp_path)

    clip_path = tmp_path / "custom_clip"
    clip_path.write_text("")
    clip_path.chmod(0o755)
    monkeypatch.setenv("IMOD_CLIP", str(clip_path))

    input_path = tmp_path / "sample.mrc"
    with mrcfile.new(input_path, overwrite=True) as handle:
        handle.set_data(np.zeros((2, 2), dtype=np.float32))

    calls = []

    def fake_runner(command, *, capture_output=False):
        calls.append((list(command), capture_output))
        return subprocess.CompletedProcess(command, 0)

    module.filter_mrc_files("0.16", command_runner=fake_runner)

    assert len(calls) == 1
    command, capture_output = calls[0]
    expected_output = "sample_filtered_by_0.16.mrc"
    assert command == [
        str(clip_path),
        "filter",
        "-l",
        "0.16",
        input_path.name,
        expected_output,
    ]
    assert capture_output is False


def test_filter_mrc_which_lookup(tmp_path, monkeypatch):
    """Verify which() lookup is used when IMOD_CLIP is not set."""
    module = _load_filter_module()
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("IMOD_CLIP", raising=False)

    input_path = tmp_path / "sample.mrc"
    with mrcfile.new(input_path, overwrite=True) as handle:
        handle.set_data(np.zeros((2, 2), dtype=np.float32))

    calls = []

    def fake_runner(command, *, capture_output=False):
        calls.append((list(command), capture_output))
        if command == ["which", "clip"]:
            return subprocess.CompletedProcess(
                command,
                0,
                stdout="/usr/local/bin/clip\n",
                stderr="",
            )
        return subprocess.CompletedProcess(command, 0)

    module.filter_mrc_files("0.33", command_runner=fake_runner)

    assert calls[0] == (["which", "clip"], True)

    command, capture_output = calls[1]
    expected_output = "sample_filtered_by_0.33.mrc"
    assert command == [
        "/usr/local/bin/clip",
        "filter",
        "-l",
        "0.33",
        input_path.name,
        expected_output,
    ]
    assert capture_output is False
