"""Tests for template_input/replace_defocus_in_star.py behaviour."""

from __future__ import annotations

import importlib.util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "utils" / "replace_defocus_in_star.py"


def _load_script_module(module_name: str):
    """Import the helper script as a module for testing."""
    spec = importlib.util.spec_from_file_location(module_name, SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)  # type: ignore[call-arg]
    return module


def test_replace_defocus_updates_expected_columns(tmp_path, monkeypatch):
    """Ensure defocus columns are replaced and STAR_DIR overrides work."""
    monkeypatch.setenv("STAR_DIR", str(tmp_path))
    module = _load_script_module("replace_defocus_in_star_template")

    # Confirm the module captured the env override at import time.
    assert module.STAR_ROOT == str(tmp_path)

    star_contents = (
        "# Example header\n"
        "#    POS something\n"
        "   1 0 0 0 0 0 10000 11000 0\n"
        "   2 0 0 0 0 0 12000 13000 0\n"
    )
    star_path = tmp_path / "tilt_example.star"
    star_path.write_text(star_contents, encoding="utf-8")

    module.replace_defocus_in_star(15000, str(tmp_path / "tilt_*.star"))

    updated_lines = star_path.read_text(encoding="utf-8").splitlines()
    # Header remains unchanged
    assert updated_lines[0] == "# Example header"
    assert updated_lines[1].startswith("#    POS")

    for line in updated_lines[2:]:
        columns = line.split()
        assert columns[6] == "15000"
        assert columns[7] == "15000"
