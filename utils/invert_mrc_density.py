import os
import subprocess
from pathlib import Path


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


def resolve_relion_binary():
    return os.environ.get("RELION_IMAGE_HANDLER", "relion_image_handler")


def invert_mrc_density(mrc_filename, *, command_runner=run_command):
    output_file_name = mrc_filename[:-4] + "_inverted.mrc"

    if os.path.exists(output_file_name):
        print(str(output_file_name) + " already exists.")
        return 1

    command = [
        resolve_relion_binary(),
        "--i",
        str(mrc_filename),
        "--o",
        str(output_file_name),
        "--multiply_constant",
        "-1",
    ]

    result = command_runner(command)
    if result.returncode != 0:
        print(f"inverting density failed for {mrc_filename}")
        return output_file_name

    return output_file_name


def invert_directory(directory: str | Path, *, command_runner=run_command) -> None:
    """Invert density for all MRC files in the given directory, renaming REC files if needed."""

    path = Path(directory)
    processed = False
    for mrc_path in path.glob("*.mrc"):
        invert_mrc_density(str(mrc_path), command_runner=command_runner)
        processed = True

    if processed:
        return

    for rec_path in path.glob("*.rec"):
        mrc_path = path / f"{rec_path.stem}_rec.mrc"
        rec_path.rename(mrc_path)
        result = invert_mrc_density(str(mrc_path), command_runner=command_runner)
        if result == 1:
            print("invert failed for", mrc_path)


if __name__ == "__main__":
    invert_directory(Path.cwd())
