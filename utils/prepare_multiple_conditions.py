"""Materialize per-PDB, per-condition input folders from a conditions file.

Environment variables
---------------------
FORCE : if set to "1" allows clobbering an existing each_input_folders directory.
GRIPTOMO_TEMPLATE_ROOT : optional override pointing to a directory containing
    template_input assets (defaults to the bundled files in ``utils``).

Parameters
----------
try_these_conditions : str
    Path to a file listing keyword followed by values (doseval, frameval, thickval, phaseval, defocusval).
folder_of_many_pdb_inputs : str
    Directory containing one or more .pdb files to expand into condition grids.

Returns
-------
bool
    True if at least one PDB was found and processed; False otherwise.
"""

import glob
import os
import shutil
import sys
import time
from pathlib import Path

MODULE_ROOT = Path(__file__).resolve().parent
DEFAULT_TEMPLATE_ENTRIES = (
    MODULE_ROOT / "template_simulate-tilt-noise.sh",
    MODULE_ROOT / "template_runscript.py",
    MODULE_ROOT / "replace_defocus_in_star.py",
    MODULE_ROOT / "run.sbatch",
    MODULE_ROOT / "angles.txt",
    MODULE_ROOT / "stars",
)

starting_dir = os.getcwd()


def _resolve_template_source(explicit: str | None) -> Path | None:
    """Resolve a template_input directory supplied via argument, env, or legacy global."""

    candidates = [
        explicit,
        os.environ.get("GRIPTOMO_TEMPLATE_ROOT"),
        globals().get("template_input_path"),
    ]
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate).expanduser()
        if not path.is_dir():
            raise FileNotFoundError(f"Template input root does not exist: {candidate}")
        return path
    return None


def _populate_template_input(target: Path, template_root: Path | None) -> None:
    """Copy template_input assets into *target* from either an override or bundled defaults."""

    if template_root is not None:
        shutil.copytree(template_root, target)
        return

    target.mkdir(parents=True, exist_ok=True)
    for entry in DEFAULT_TEMPLATE_ENTRIES:
        if not entry.exists():
            continue
        destination = target / entry.name
        if entry.is_dir():
            shutil.copytree(entry, destination)
        else:
            shutil.copy2(entry, destination)


def prepare_many_inputs(
    try_these_conditions,
    folder_of_many_pdb_inputs,
    template_input_root: str | None = None,
):
    output_folder = "each_input_folders"
    if os.path.exists(output_folder):
        if os.environ.get("FORCE") == "1":
            shutil.rmtree(output_folder)
        else:
            raise FileExistsError(
                f"{output_folder} already exists; set FORCE=1 to overwrite existing inputs"
            )

    os.makedirs(output_folder, exist_ok=True)

    pdb_exists = False
    glob_these = os.path.join(folder_of_many_pdb_inputs, "*.pdb")
    template_root = _resolve_template_source(template_input_root)

    for each_input_pdb in glob.glob(glob_these):
        pdb_exists = True
        print(f"processing pdb: {each_input_pdb}")

        each_input_foldername = os.path.basename(each_input_pdb[:-4])
        new_template_input_path = (
            Path(output_folder) / each_input_foldername / "template_input"
        )
        new_template_input_path.parent.mkdir(parents=True, exist_ok=True)
        _populate_template_input(new_template_input_path, template_root)
        shutil.copy(each_input_pdb, new_template_input_path)

        template_sh = new_template_input_path / "template_simulate-tilt-noise.sh"
        output_flag = new_template_input_path / "simulate-tilt-noise.sh"
        with template_sh.open() as f_in, output_flag.open("w") as f_out:
            for line in f_in:
                new_line = line.replace(
                    "replace_here", os.path.join("..", os.path.basename(each_input_pdb))
                )
                f_out.write(new_line)

        # Build value lists per keyword.
        list_doseval: list[str] = []
        list_frameval: list[str] = []
        list_thickval: list[str] = []
        list_phaseval: list[str] = []
        list_defocusval: list[str] = []
        with open(try_these_conditions) as f_cond_in:
            for line in f_cond_in:
                line_list = line.split()
                if not line_list:
                    continue
                key, *values = line_list
                target = None
                if key == "doseval":
                    target = list_doseval
                elif key == "frameval":
                    target = list_frameval
                elif key == "thickval":
                    target = list_thickval
                elif key == "phaseval":
                    target = list_phaseval
                elif key == "defocusval":
                    target = list_defocusval
                if target is not None:
                    target.extend(values)

        # Cartesian expansion of all lists to build condition folder names.
        for i in range(len(list_doseval)):
            for j in range(len(list_frameval)):
                for k in range(len(list_thickval)):
                    for l_idx in range(len(list_phaseval)):
                        for m in range(len(list_defocusval)):
                            each_condition_foldername = f"{list_doseval[i]}x{list_frameval[j]}e_{list_thickval[k]}t_{list_phaseval[l_idx]}p_{list_defocusval[m]}d"
                            final_input_path = (
                                Path(output_folder)
                                / each_input_foldername
                                / each_condition_foldername
                            )

                            shutil.copytree(new_template_input_path, final_input_path)

                            template_runscript = (
                                final_input_path / "template_runscript.py"
                            )
                            final_runscript = final_input_path / "runscript.py"
                            with (
                                template_runscript.open() as f_in,
                                final_runscript.open("w") as f_out,
                            ):
                                for line in f_in:
                                    split_line = line.split()
                                    if len(split_line) < 2:
                                        f_out.write(line)
                                        continue
                                    value_map = {
                                        "doseval": list_doseval[i],
                                        "frameval": list_frameval[j],
                                        "thickval": list_thickval[k],
                                        "phaseval": list_phaseval[l_idx],
                                        "defocusval": list_defocusval[m],
                                    }
                                    key = split_line[1]
                                    replacement = value_map.get(key)
                                    if replacement:
                                        f_out.write(
                                            line.replace("replace_this", replacement)
                                        )
                                    else:
                                        f_out.write(line)

    return pdb_exists


if __name__ == "__main__":
    ########## <begin> basic checks
    if not hasattr(sys, "version_info") or sys.version_info < (3, 7):
        raise ValueError("Script requires Python 3.7 or higher!")
    args = sys.argv[:]
    if len(args) < 3:
        print(f"len(args):{len(args)}")
        print(
            "Specify                               <try_these_conditions>         <folder_of_many_pdb_inputs>"
        )
        print("For example,\n")
        print(
            "python prepare_multiple_conditions.py input/try_these_conditions.txt input/many_pdb_inputs"
        )
        exit(1)
    ########## <end> parameter checks

    py_code = args[0]
    try_these_conditions = os.path.abspath(args[1])
    folder_of_many_pdb_inputs = os.path.abspath(args[2])

    start_time = time.time()
    pdb_exists = prepare_many_inputs(
        try_these_conditions,
        folder_of_many_pdb_inputs,
        os.environ.get("GRIPTOMO_TEMPLATE_ROOT"),
    )
    if not pdb_exists:
        print(f"No pdb files found in {folder_of_many_pdb_inputs}")
        exit(1)
    end_time = time.time()

    elapsed_time = end_time - start_time
    if elapsed_time < 3600:
        print(f"Preparing many inputs elapsed: {elapsed_time / 60:.1f} minutes")
    else:
        print(f"Preparing many inputs elapsed: {elapsed_time / 3600:.1f} hours")
