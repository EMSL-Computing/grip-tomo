#!/usr/bin/env python
"""Merge simulated MRC files into a single MRCS stack."""

from __future__ import annotations

import os
import shlex
import shutil
import time
from pathlib import Path

import mrcfile
import numpy as np
import pandas as pd


def check_which_mrc_file_can_be_merged_by_averaging(
    df_of_all_mrc_files_sorted_by_angle: pd.DataFrame,
) -> pd.DataFrame:
    """Filter out MRC files that are empty or exhibit zero-mean density."""
    records: list[dict[str, float | str]] = []
    for _, row in df_of_all_mrc_files_sorted_by_angle.iterrows():
        mrc_filename = row["mrc_filename"]
        file_size = os.path.getsize(mrc_filename)
        if file_size == 0:
            raise SystemExit(f"{mrc_filename} has zero file size; aborting merge.")

        with mrcfile.open(mrc_filename) as each_mrc:
            if float(np.mean(each_mrc.data)) == 0:
                print(
                    "warning: mean density is zero for "
                    f"{mrc_filename}; EMAN2 may skip this tilt during merge."
                )
            else:
                records.append(
                    {
                        "mrc_filename": str(mrc_filename),
                        "angle": float(row["angle"]),
                    }
                )

    df_of_proper_mrc_files = pd.DataFrame(records, columns=["mrc_filename", "angle"])
    if not df_of_proper_mrc_files.empty:
        df_of_proper_mrc_files.sort_values(by="angle", ascending=True, inplace=True)
        df_of_proper_mrc_files.reset_index(inplace=True, drop=True)
    df_of_proper_mrc_files.to_csv(
        "df_of_proper_mrc_files_sorted_by_angle.csv", index=True
    )
    return df_of_proper_mrc_files


def when_all_mrc_files_are_proper_to_be_merged(
    df_of_all_mrc_files_sorted_by_angle: pd.DataFrame, e2proc2d_bin: str
) -> None:
    """Invoke EMAN2's e2proc2d to merge validated MRC files into an MRCS stack."""
    command_parts = [shlex.quote(e2proc2d_bin)]
    for mrc_filename in df_of_all_mrc_files_sorted_by_angle["mrc_filename"]:
        file_size = os.path.getsize(mrc_filename)
        if file_size == 0:
            raise SystemExit(f"{mrc_filename} has zero file size; aborting merge.")
        command_parts.append(shlex.quote(str(mrc_filename)))

    base_mrc_filename = os.path.basename(mrc_filename)

    n = base_mrc_filename.count("_")
    up_to_here = -1  # starting
    for i in range(0, n):
        up_to_here = base_mrc_filename.find("_", up_to_here + 1)

    final_mrcs_filename = base_mrc_filename[:up_to_here] + ".mrcs"
    command_parts.append(shlex.quote(str(final_mrcs_filename)))
    command = " ".join(command_parts)
    print("mrc files merging command")
    print(command)
    os.system(command)
    # EMAN2 skips frames that trigger the "sigma = 0" warning, so inputs should
    # be screened before reaching this point.


def _locate_e2proc2d() -> str:
    """Resolve the e2proc2d binary using env overrides or PATH lookup."""
    override = os.environ.get("EMAN2_E2PROC2D")
    if override:
        override_path = Path(override)
        if override_path.exists():
            return str(override_path)
        resolved_override = shutil.which(override)
        if resolved_override:
            return resolved_override
        raise SystemExit(
            "EMAN2_E2PROC2D is set but the specified executable was not found: "
            f"{override}"
        )

    resolved = shutil.which("e2proc2d.py")
    if resolved:
        return resolved

    raise SystemExit(
        "e2proc2d.py is not available. Set EMAN2_E2PROC2D or activate the EMAN2 environment."
    )


if __name__ == "__main__":
    start_time = time.time()

    ## <begin> basic check of requirement
    e2proc2d_bin = _locate_e2proc2d()
    print(f"using e2proc2d executable: {e2proc2d_bin}")
    ## <end> basic check of requirement

    df_of_all_mrc_files_sorted_by_angle_filename = (
        "df_of_all_mrc_files_sorted_by_angle.csv"
    )

    colnames = ["mrc_filename", "angle"]
    df_of_all_mrc_files_sorted_by_angle = pd.read_csv(
        df_of_all_mrc_files_sorted_by_angle_filename, names=colnames, skiprows=1
    )

    # df_of_proper_mrc_files = check_which_mrc_file_can_be_merged_by_trying(df_of_all_mrc_files_sorted_by_angle)
    df_of_proper_mrc_files = check_which_mrc_file_can_be_merged_by_averaging(
        df_of_all_mrc_files_sorted_by_angle
    )

    if df_of_proper_mrc_files.empty:
        raise SystemExit("No suitable MRC files were detected; aborting merge.")

    when_all_mrc_files_are_proper_to_be_merged(df_of_proper_mrc_files, e2proc2d_bin)

    duration = time.time() - start_time
    print(f"reconstruction took {duration:.2f}s")
