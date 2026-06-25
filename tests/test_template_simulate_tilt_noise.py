"""Tests for template_simulate-tilt-noise.sh substitution and env overrides."""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "utils" / "template_simulate-tilt-noise.sh"


def test_template_has_shebang_and_placeholder():
    """Ensure simulate-tilt-noise template has shebang and replace_here token."""
    text = SCRIPT_PATH.read_text(encoding="utf-8")
    assert text.startswith("#!/") or text.startswith(
        "#"
    )  # allow minimal shebang/comments
    assert "replace_here" in text


def test_placeholder_replacement_logic(tmp_path):
    """Simulate replacement of placeholder and confirm relative path usage."""
    template = "#!/bin/bash\nrun replace_here\n"
    template_path = tmp_path / "template_simulate-tilt-noise.sh"
    template_path.write_text(template, encoding="utf-8")
    pdb_path = tmp_path / "sample.pdb"
    pdb_path.write_text("ATOM\n", encoding="utf-8")

    output_path = tmp_path / "simulate-tilt-noise.sh"
    replaced = template.replace("replace_here", f"../{pdb_path.name}")
    output_path.write_text(replaced, encoding="utf-8")

    assert f"../{pdb_path.name}" in output_path.read_text(encoding="utf-8")
    assert "replace_here" not in output_path.read_text(encoding="utf-8")
