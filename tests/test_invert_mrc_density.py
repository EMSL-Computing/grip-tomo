import subprocess
import importlib


def _load_invert_module():
    return importlib.import_module("utils.invert_mrc_density")


def test_invert_mrc_density_uses_env_override(tmp_path, monkeypatch):
    """Ensure RELION env override is honored when inverting MRC density."""
    module = _load_invert_module()
    monkeypatch.chdir(tmp_path)

    relion_path = tmp_path / "relion_image_handler_custom"
    relion_path.write_text("")
    monkeypatch.setenv("RELION_IMAGE_HANDLER", str(relion_path))

    input_path = tmp_path / "input.mrc"
    input_path.write_bytes(b"")

    calls = []

    def fake_runner(command, *, capture_output=False):
        calls.append((list(command), capture_output))
        return subprocess.CompletedProcess(command, 0)

    output_path = module.invert_mrc_density(str(input_path), command_runner=fake_runner)

    assert len(calls) == 1
    command, capture_output = calls[0]
    expected_output = str(tmp_path / "input_inverted.mrc")
    assert command == [
        str(relion_path),
        "--i",
        str(input_path),
        "--o",
        expected_output,
        "--multiply_constant",
        "-1",
    ]
    assert capture_output is False
    assert output_path == expected_output


def test_invert_mrc_density_skips_when_output_exists(tmp_path):
    """Skip inversion when the output inverted MRC already exists."""
    module = _load_invert_module()

    input_path = tmp_path / "input.mrc"
    input_path.write_bytes(b"")
    output_path = tmp_path / "input_inverted.mrc"
    output_path.write_bytes(b"")

    calls = []

    def fake_runner(command, *, capture_output=False):
        calls.append((list(command), capture_output))
        return subprocess.CompletedProcess(command, 0)

    result = module.invert_mrc_density(str(input_path), command_runner=fake_runner)

    assert result == 1
    assert calls == []


def test_invert_mrc_density_logs_failure(tmp_path, monkeypatch, capsys):
    """Print helpful hint and return expected filename on relion failure."""
    module = _load_invert_module()
    input_path = tmp_path / "input.mrc"
    input_path.write_bytes(b"")

    def fake_run(command, *, check=False, capture_output=False, text=False):
        return subprocess.CompletedProcess(command, 2)

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    output_path = module.invert_mrc_density(str(input_path))

    captured = capsys.readouterr()
    assert "return code: 2" in captured.out
    assert output_path == str(tmp_path / "input_inverted.mrc")
