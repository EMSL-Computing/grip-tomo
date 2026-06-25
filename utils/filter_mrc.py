import glob
import os
import platform
import subprocess
import sys
from pathlib import Path

import mrcfile


def run_command(command, *, capture_output=False):
    printable = " ".join(str(part) for part in command)
    print(printable)
    completed = subprocess.run(
        [str(part) for part in command],
        check=False,
        capture_output=capture_output,
        text=True if capture_output else False,
    )
    print(f"return code: {completed.returncode}")
    return completed


def resolve_clip_binary(command_runner=run_command):
    env_value = os.environ.get("IMOD_CLIP")
    if env_value:
        clip_path = Path(env_value).expanduser()
        if clip_path.is_file() and os.access(clip_path, os.X_OK):
            print(str(clip_path))
            return str(clip_path)
        result = command_runner(["which", str(env_value)], capture_output=True)
        if result.returncode == 0 and result.stdout.strip():
            clip_candidate = result.stdout.strip()
            print(clip_candidate)
            return clip_candidate
        print("clip is not found")
        sys.exit(1)
    result = command_runner(["which", "clip"], capture_output=True)
    if result.returncode != 0 or not result.stdout.strip():
        print("clip is not found")
        sys.exit(1)
    clip_path = result.stdout.strip()
    print(clip_path)
    return clip_path


def filter_mrc_files(low_pass_level, *, command_runner=run_command):
    clip_binary = resolve_clip_binary(command_runner)
    patterns = ("*.mrc", "*.rec")
    files = []
    for pattern in patterns:
        files.extend(glob.glob(pattern))

    for input_mrc in files:
        with mrcfile.open(input_mrc):
            output_mrc = (
                input_mrc[:-4] + "_filtered_by_" + str(low_pass_level) + input_mrc[-4:]
            )
        command_runner(
            [
                clip_binary,
                "filter",
                "-l",
                str(low_pass_level),
                str(input_mrc),
                str(output_mrc),
            ]
        )


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    instructions = """
0.05  low pass level -> remove data better than 20 Angstrom resolution

0.10  low pass level -> remove data better than 10 Angstrom resolution
0.125 low pass level -> remove data better than 8 Angstrom resolution
0.16  low pass level -> remove data better than 6 Angstrom resolution
0.20  low pass level -> remove data better than 5 Angstrom resolution
0.25  low pass level -> remove data better than 4 Angstrom resolution
0.33  low pass level -> remove data better than 3 Angstrom resolution
0.4   low pass level -> remove data better than 2.5 Angstrom resolution"""
    print(instructions)

    if not argv:
        print("Specify                                    <low pass level>")
        print("For example,\n")
        print("python After_cistem_till_reconstruction.py 0.16")
        return 1

    low_pass_level = argv[0]
    filter_mrc_files(low_pass_level)
    return 0


if __name__ == "__main__":
    sys.exit(main())
