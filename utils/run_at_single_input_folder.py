#!/usr/bin/env python
# coding: utf-8

# run @.../input
# assumed that a user ran cistem at input folder

import os
from pathlib import Path

_output_env = os.environ.get("GRIPTOMO_OUTPUT_DIR")
OUTPUT_DIR = Path(_output_env) if _output_env else Path("../output")


def move_outputs(workdir: Path, output_dir: Path = OUTPUT_DIR) -> None:
    """Move cisTEM and simulate_* outputs from workdir to output_dir, preserving other files.

    Parameters
    ----------
    workdir : Path
        Directory containing cisTEM run outputs.
    output_dir : str
        Target directory (may be relative) to receive selected outputs.
    """
    target = output_dir
    if not target.is_absolute():
        if str(target).startswith("../"):
            # place outside workdir (sibling to workdir)
            target = (workdir.parent / target.parts[-1]).resolve()
        elif ".." not in str(target):
            target = (workdir / target).resolve()
    target.mkdir(parents=True, exist_ok=True)

    patterns = ["cistem.", "simulate_"]
    for item in workdir.iterdir():
        name = item.name
        if any(name.startswith(pat) for pat in patterns):
            dest = target / name
            item.replace(dest)


if __name__ == "__main__":
    starting_dir = Path.cwd()
    print(f"starting_dir:{starting_dir}")
    move_outputs(starting_dir, OUTPUT_DIR)
    print("See", OUTPUT_DIR)
