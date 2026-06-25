"""Smoke tests for template_input/template_runscript.py integrity."""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "utils" / "template_runscript.py"


def _load_template_text() -> str:
    """Read the template runscript contents for inspection."""
    return SCRIPT_PATH.read_text(encoding="utf-8")


def _expected_angles() -> list[float]:
    """Return the canonical list of tilt angles used in the template."""
    values = []
    angle = -90.0
    while angle <= 90.0 + 1e-9:
        values.append(round(angle, 1))
        angle += 0.5
    return values


def test_template_contains_placeholder_loops():
    """Ensure template still exposes replace_this placeholders for grid search."""
    text = _load_template_text()
    assert "for doseval in [replace_this]" in text
    assert "for frameval in [replace_this]" in text
    assert "for thickval in [replace_this]" in text
    assert "for defocusval in [replace_this]" in text


def test_template_angles_sequence_matches_expected(tmp_path):
    """Ensure template provides the canonical 0.5 step tilt angle sequence."""
    text = _load_template_text()
    expected = _expected_angles()

    angles: list[float] = []
    marker = "for angleval in ["
    if marker in text:
        segment = text.split(marker, 1)[1].split(":", 1)[0]
        segment = segment.split("]", 1)[0]
        for token in segment.split(","):
            stripped = token.strip()
            if stripped:
                angles.append(float(stripped))
    else:
        angles_path = SCRIPT_PATH.with_name("angles.txt")
        contents = angles_path.read_text(encoding="utf-8")
        for line in contents.splitlines():
            if line.strip():
                angles.append(float(line.strip()))

    assert angles == expected


def test_template_invokes_expected_helpers():
    """Ensure template points to helper scripts with relative paths."""
    text = _load_template_text()
    assert "../replace_defocus_in_star.py" in text
    assert "../simulate-tilt-noise.sh" in text
