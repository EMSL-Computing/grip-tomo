#!/usr/bin/env python
"""Parallel dispatcher for running single-folder reconstruction across condition directories.

Environment variables
---------------------
GRIPTOMO_PYTHON_BIN : Override Python interpreter used to invoke the single-folder script.
GRIPTOMO_MAX_WORKERS : Optional integer limiting ProcessPoolExecutor workers.
GRIPTOMO_SINGLE_FOLDER_SCRIPT : Optional override for at_single_folder.py location.
"""

import concurrent.futures
import os
import subprocess
import sys
import time
from pathlib import Path

DEFAULT_PYTHON_BIN = sys.executable


def process_directory(
    d_simulate___, run_this_py, min_tilt, max_tilt, increment, low_pass_level
):
    """Invoke the single-folder reconstruction script for one simulate_* directory."""
    python_bin = os.environ.get("GRIPTOMO_PYTHON_BIN", DEFAULT_PYTHON_BIN)
    command = [
        python_bin,
        str(run_this_py),
        str(min_tilt),
        str(max_tilt),
        str(increment),
        str(low_pass_level),
    ]
    print(f"Current working subfolder: {d_simulate___}")
    completed = subprocess.run(command, check=False, cwd=d_simulate___)
    if completed.returncode != 0:
        print(f"Command exited with code {completed.returncode} in {d_simulate___}")


### end of def process_directory(d_simulate___, run_this_py, min_tilt, max_tilt, increment, low_pass_level):


def modified_run_at_many_folders(
    run_this_py, min_tilt, max_tilt, increment, low_pass_level
):
    """Traverse each_input_folders and submit reconstruction tasks for each simulate_* directory."""
    futures = []
    max_workers_env = os.environ.get("GRIPTOMO_MAX_WORKERS")
    executor_kwargs = {}
    if max_workers_env:
        try:
            executor_kwargs["max_workers"] = int(max_workers_env)
        except ValueError:
            print(
                f"Invalid GRIPTOMO_MAX_WORKERS value: {max_workers_env}; using default"
            )
    with concurrent.futures.ProcessPoolExecutor(**executor_kwargs) as executor:
        for file in os.listdir(starting_dir):
            each_input_folder = os.path.join(starting_dir, file)
            if os.path.isdir(each_input_folder):
                print(f"(level 1) each_input_folder:{each_input_folder}")
                for each_condition_folder in os.listdir(each_input_folder):
                    print(f"(level 2) each_condition_folder:{each_condition_folder}")

                    if each_condition_folder == "template_input":
                        continue

                    output_dir = os.path.join(
                        each_input_folder, each_condition_folder, "output"
                    )
                    print(f"output_dir:{output_dir}")
                    if not os.path.isdir(output_dir):
                        print("output folder is not found.")
                        print(
                            "Run griptomoml/scripts/prepare_data/mrc/1_prepare_2d/w_cistem/1_pdb_to_2d_images/multiple_pdbs/either_multiple_pdbs_or_1_pdb_input/3_process_output first"
                        )
                        sys.exit(1)
                    os.chdir(output_dir)
                    for file in os.listdir(output_dir):
                        d_simulate___ = os.path.join(output_dir, file)
                        if os.path.isdir(d_simulate___):
                            future = executor.submit(
                                process_directory,
                                d_simulate___,
                                run_this_py,
                                min_tilt,
                                max_tilt,
                                increment,
                                low_pass_level,
                            )
                            futures.append(future)

        # Wait for all tasks to complete
        for future in concurrent.futures.as_completed(futures):
            future.result()


### end of def modified_run_at_many_folders(run_this_py, min_tilt, max_tilt, increment, low_pass_level):


if __name__ == "__main__":
    # Version guard & dependency presence checks.
    if not hasattr(sys, "version_info") or sys.version_info < (3, 7):
        raise ValueError("Script requires Python 3.7 or higher!")

    try:
        subprocess.check_output("which e2pdb2mrc.py", shell=True)
    except Exception:
        sys.exit(
            "e2pdb2mrc.py not found in PATH. Activate EMAN2 environment or add its bin directory (override via PATH)."
        )

    try:
        subprocess.check_output("which header", shell=True)
    except Exception:
        sys.exit(
            "IMOD 'header' executable not found; source IMOD environment so header is on PATH."
        )

    try:
        subprocess.check_output("which relion_image_handler", shell=True)
    except Exception:
        sys.exit(
            "relion_image_handler not found; source RELION environment so relion_image_handler is on PATH."
        )

    print_this = """    0.10  low pass level -> remove above 10 Angstrom\n\
    0.125 low pass level -> remove above 8 Angstrom\n\
    0.16  low pass level -> remove above 6 Angstrom\n\
    0.20  low pass level -> remove above 5 Angstrom\n\
    0.25  low pass level -> remove above 4 Angstrom\n\
    0.33  low pass level -> remove above 3 Angstrom\n\
    0.4   low pass level -> remove above 2.5 Angstrom"""
    print(print_this)

    args = sys.argv[:]

    if len(args) < 5:
        print(
            "Specify                                    <min_tilt in integer> <max_tilt in integer> <increment in integer> <low pass level>"
        )
        print("For example,")
        print(
            "python After_cistem_till_reconstruction.py -60                   60                    3                      0.16"
        )
    sys.exit(1)
    ########## <end> basic checks
    starting_dir = os.getcwd()

    py_code = args[0]

    code_location = os.path.dirname(os.path.abspath(py_code))
    run_this_py = os.environ.get("GRIPTOMO_SINGLE_FOLDER_SCRIPT")
    if run_this_py:
        run_this_py = os.path.abspath(run_this_py)
    else:
        run_this_py = str(Path(__file__).resolve().with_name("at_single_folder.py"))

    min_tilt = int(args[1])
    max_tilt = int(args[2])
    increment = float(args[3])
    low_pass_level = float(args[4])

    start_time = time.time()
    modified_run_at_many_folders(
        run_this_py, min_tilt, max_tilt, increment, low_pass_level
    )

    print(
        "2D->3D for all input pdb files and conditions took --- %s minutes ---"
        % round((time.time() - start_time) / 60, 2)
    )
