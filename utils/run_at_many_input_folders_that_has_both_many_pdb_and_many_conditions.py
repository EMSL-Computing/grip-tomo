#!/usr/bin/env python
"""Organize cisTEM outputs into per-condition input/output directories.

Environment variables
---------------------
GRIPTOMO_SINGLE_PROCESSOR : Optional path override for the single-folder processor script.
"""

import glob
import os
import shutil
import sys

OUTPUT_MOVE_PATTERNS = (
    "cisTEM.*",
    "simulate_*x*e_*t_*p",
)


def run_at_many_folders(run_this_py):
    """Walk each_input_folders and move cisTEM outputs into structured folders.

    Parameters
    ----------
    run_this_py : str
        Path to helper script (currently unused but retained for compatibility).
    """
    for folder_candidate in os.listdir(starting_dir):
        each_pdb_input_folder = os.path.join(starting_dir, folder_candidate)
        if os.path.isdir(each_pdb_input_folder):
            for each_condition_folder in os.listdir(each_pdb_input_folder):
                if "template_input" not in each_condition_folder:
                    os.chdir(os.path.join(each_pdb_input_folder, each_condition_folder))
                    if not os.path.isdir("output"):
                        if (
                            len(
                                [
                                    each_file
                                    for each_file in os.listdir(os.getcwd())
                                    if each_file.startswith("cisTEM")
                                ]
                            )
                            > 0
                        ):
                            os.makedirs("output", exist_ok=True)
                            for pattern in OUTPUT_MOVE_PATTERNS:
                                for match in glob.glob(pattern):
                                    destination = os.path.join(
                                        "output", os.path.basename(match)
                                    )
                                    shutil.move(match, destination)

                            if not os.path.isdir("input"):
                                os.mkdir("input")
                                for each_file in os.listdir(os.getcwd()):
                                    if each_file != "input":
                                        shutil.move(
                                            each_file, os.path.join("input", each_file)
                                        )

                            moved_output_path = os.path.join("input", "output")
                            if os.path.exists(moved_output_path):
                                shutil.move(moved_output_path, ".")

                    else:
                        print(f"output folder already exists: {os.getcwd()}")
                        continue
    print("See output folder at each condition and input pdb file")


###### end of def run_at_many_folders():

if __name__ == "__main__":
    # Basic version guard.
    if not hasattr(sys, "version_info") or sys.version_info < (3, 7):
        raise ValueError("Script requires Python 3.7 or higher!")
    args = sys.argv[:]
    # Determine starting directory and resolve helper script path.
    starting_dir = os.getcwd()

    py_code = args[0]

    code_location = os.path.dirname(os.path.abspath(py_code))
    run_this_py = os.environ.get("GRIPTOMO_SINGLE_PROCESSOR")
    if run_this_py:
        run_this_py = os.path.abspath(run_this_py)
    # else: legacy fallback removed; require explicit env if needed.

    run_at_many_folders(run_this_py)
