"""Smoke tests for utils.sliding_window."""

from __future__ import annotations

import csv
from pathlib import Path

import mrcfile
import numpy as np

from utils import sliding_window


def _write_mrc(path: Path, data: np.ndarray) -> None:
    with mrcfile.new(str(path), overwrite=True) as mrc:
        mrc.set_data(data.astype(np.float32))
        mrc.header.origin = (0.0, 0.0, 0.0)
        mrc.voxel_size = (1.0, 1.0, 1.0)


def test_generate_subvolumes_smoke(tmp_path, monkeypatch):
    """Ensure generate_subvolumes produces expected outputs and metadata."""

    volume = np.arange(27, dtype=np.float32).reshape(3, 3, 3)
    input_path = tmp_path / "input.mrc"
    _write_mrc(input_path, volume)

    monkeypatch.chdir(tmp_path)

    sub_size = (2, 2, 2)
    step = (1, 1, 1)
    outputs = sliding_window.generate_subvolumes(
        str(input_path),
        "sub",
        subvolume_size=sub_size,
        step_size=step,
    )

    assert len(outputs) == 8
    for output_name in outputs:
        out_path = tmp_path / output_name
        assert out_path.exists(), f"missing subvolume {output_name}"
        with mrcfile.open(out_path) as sub_mrc:
            assert sub_mrc.data.shape == sub_size

    csv_log = tmp_path / f"{sub_size}_{step}.csv"
    assert csv_log.exists()
    with csv_log.open(newline="") as handle:
        reader = csv.reader(handle)
        rows = list(reader)
    assert len(rows) == 1 + len(outputs)
    header = rows[0]
    assert header == [
        "file name",
        "x-start",
        "x-end",
        "y-start",
        "y-end",
        "z-start",
        "z-end",
    ]
