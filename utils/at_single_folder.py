#!/usr/bin/env python
"""Run the post-cisTEM reconstruction pipeline in a single simulation folder."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd


@dataclass(frozen=True)
class PipelinePaths:
    """Container for helper script locations."""

    scripts_root: Path
    merge_script: Path
    phase_template: Path
    reconstruct_script: Path
    # filter & invert handled via direct module imports (utils.filter_mrc / utils.invert_mrc_density)


def generate_angles(min_tilt: float, max_tilt: float, increment: float) -> list[float]:
    """Return the ordered list of tilt angles for the reconstruction range."""

    angles: list[float] = []
    current = float(min_tilt)
    epsilon = max(increment * 0.01, 1e-6)
    while current <= max_tilt + epsilon:
        angles.append(round(current, 6))
        current += increment
    return angles


def shlex_join(cmd: Sequence[str]) -> str:
    """Join a command list into a printable shell string."""

    try:
        return shlex.join(cmd)
    except AttributeError:  # Python <3.8 compatibility
        return " ".join(shlex.quote(part) for part in cmd)


def run_command(
    cmd: Sequence[str], cwd: Path | None = None
) -> subprocess.CompletedProcess:
    """Execute a command, printing output and raising on failure."""

    display = shlex_join(list(map(str, cmd)))
    if cwd:
        print(f"[cwd={cwd}] {display}")
    else:
        print(display)

    result = subprocess.run(
        list(map(str, cmd)),
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
    )

    if result.stdout:
        print(result.stdout.rstrip())
    if result.stderr:
        print(result.stderr.rstrip(), file=sys.stderr)

    if result.returncode != 0:
        raise SystemExit(
            f"command failed with exit code {result.returncode}: {display}"
        )

    return result


def mrc_to_mrcs(
    start_dir: Path,
    angles: Iterable[float],
    paths: PipelinePaths,
    python_bin: str,
) -> None:
    """Prepare the sorted MRC listing and merge into an MRCS stack."""

    for existing_mrcs in start_dir.glob("*.mrcs"):
        existing_mrcs.unlink()

    records: list[dict[str, float | str]] = []
    for mrc_file in sorted(start_dir.glob("*.mrc")):
        if "noise" not in mrc_file.name:
            continue

        token = mrc_file.name.split("_")[-1]
        angle_str = token.replace("deg.mrc", "")
        try:
            angle = float(angle_str)
        except ValueError:
            continue

        records.append({"mrc_filename": str(mrc_file), "angle": round(angle, 6)})

    if not records:
        raise SystemExit(
            "No candidate noise MRC files found; run this script inside a simulation folder."
        )

    df = pd.DataFrame(records, columns=["mrc_filename", "angle"])
    valid_angles = {round(value, 6) for value in angles}
    df = df[df["angle"].isin(valid_angles)]
    df.sort_values("angle", inplace=True)
    df.reset_index(drop=True, inplace=True)

    df_path = start_dir / "df_of_all_mrc_files_sorted_by_angle.csv"
    df.to_csv(df_path, index=True)

    run_command([python_bin, str(paths.merge_script)], cwd=start_dir)


def _write_angles_file(start_dir: Path, angles: Iterable[float]) -> Path:
    angles_path = start_dir / "angles.rawtlt"
    angles_path.unlink(missing_ok=True)
    with angles_path.open("w", encoding="utf-8") as handle:
        for value in angles:
            handle.write(f"{value}\n")
    return angles_path


def _build_defocus_file(
    base_name: str,
    start_dir: Path,
    phase_template: Path,
) -> Path:
    defocus_path = start_dir / "in_nm.defocus"
    if not phase_template.exists():
        raise SystemExit(f"Phase template not found: {phase_template}")

    defocus_value = "0"
    if "noise" in base_name:
        token = base_name.split("_")[-1]
        defocus = token.split("d")[0]
        try:
            defocus_value = str(int(int(defocus) / 10))
        except ValueError:
            defocus_value = "0"

    with (
        phase_template.open("r", encoding="utf-8") as source,
        defocus_path.open("w", encoding="utf-8") as target,
    ):
        for line in source:
            if line.startswith("#"):
                continue
            if line.startswith("12"):
                target.write(line)
                continue

            parts = line.split()
            if len(parts) < 6:
                target.write(line)
                continue

            parts[4] = defocus_value
            target.write(" ".join(parts) + "\n")

    return defocus_path


def phase_flip(
    start_dir: Path,
    angles: Iterable[float],
    paths: PipelinePaths,
    ctfphaseflip_bin: str,
) -> None:
    """Apply CTF phase flipping to each merged tilt stack."""

    angles_path = _write_angles_file(start_dir, angles)

    for mrcs_file in sorted(start_dir.glob("*.mrcs")):
        if "ctfcorrected" in mrcs_file.name:
            continue

        base_name = mrcs_file.name
        print(f"ctf phase flip {base_name}")

        phase_token = base_name.split("_")[-2]
        deg_phase = float(phase_token.split("p")[0]) * 180.0
        output_file = (
            start_dir / f"{mrcs_file.stem}_degPhase_{deg_phase}_ctfcorrected.mrcs"
        )

        defocus_file = _build_defocus_file(base_name, start_dir, paths.phase_template)

        cmd = [
            ctfphaseflip_bin,
            "-input",
            str(mrcs_file),
            "-output",
            str(output_file),
            "-angleFn",
            str(angles_path),
            "-pixelSize",
            "0.1",
            "-defTol",
            "200",
            "-volt",
            "300",
            "-cs",
            "2.7",
            "-ampContrast",
            "0.07",
            "-iWidth",
            "15",
            "-scale",
            "0.25",
            "-degPhase",
            str(deg_phase),
            "-defFn",
            str(defocus_file),
        ]

        run_command(cmd, cwd=start_dir)


def reconstruct_from_2d(start_dir: Path, paths: PipelinePaths) -> None:
    """Invoke the IMOD back-projection pipeline and clean intermediate outputs."""

    run_command(["bash", str(paths.reconstruct_script)], cwd=start_dir)

    for rec_file in start_dir.glob("*.rec"):
        if "ctfcorrected" not in rec_file.name:
            rec_file.unlink()


def apply_filters(
    rec_dir: Path, low_pass_level: float, *, filter_module=None, invert_module=None
) -> None:
    """Run filtering and inversion using the flattened utils modules."""

    if filter_module is None or invert_module is None:
        try:
            import utils.filter_mrc as filter_module  # type: ignore
            import utils.invert_mrc_density as invert_module  # type: ignore
        except ImportError as exc:  # pragma: no cover
            raise SystemExit(f"Failed to import filtering helpers: {exc}") from exc

    original_cwd = Path.cwd()
    try:
        os.chdir(rec_dir)
        filter_module.filter_mrc_files(str(low_pass_level))
        invert_module.invert_directory(Path.cwd())
    finally:
        os.chdir(original_cwd)


def restore_previous_run(start_dir: Path) -> None:
    """Move archived MRC files from a previous run back into the working folder."""

    mrc_dir = start_dir / "mrc"
    if not mrc_dir.is_dir():
        return
    for mrc_file in mrc_dir.glob("*.mrc"):
        shutil.move(str(mrc_file), start_dir / mrc_file.name)


def move_outputs(start_dir: Path, rec_dir: Path) -> None:
    """Collect generated artefacts into the reconstruction directory."""

    for pattern in ("*.mrcs", "*.csv", "*.defocus"):
        for artifact in start_dir.glob(pattern):
            shutil.move(str(artifact), rec_dir / artifact.name)

    for rec_file in start_dir.glob("*.rec"):
        renamed = rec_file.with_name(f"{rec_file.stem}_{rec_dir.name}.rec")
        rec_file.rename(renamed)
        shutil.move(str(renamed), rec_dir / renamed.name)

    for raw_file in start_dir.glob("*.rawtlt"):
        shutil.move(str(raw_file), rec_dir / raw_file.name)


def archive_remaining(start_dir: Path) -> None:
    """Move residual MRC and STAR files into their archival directories."""

    mrc_dir = start_dir / "mrc"
    mrc_dir.mkdir(exist_ok=True)
    for mrc_file in start_dir.glob("*.mrc"):
        shutil.move(str(mrc_file), mrc_dir / mrc_file.name)

    star_dir = start_dir / "star"
    star_dir.mkdir(exist_ok=True)
    for star_file in start_dir.glob("*.star"):
        shutil.move(str(star_file), star_dir / star_file.name)


def _find_util_root(script_file: Path) -> Path:
    """Locate the nearest ancestor directory named 'utils'."""

    for parent in script_file.parents:
        if parent.name == "utils":
            return parent
    raise SystemExit("Unable to locate the 'utils' directory relative to this file.")


def _resolve_path(env_var: str, default: Path) -> Path:
    candidate = os.environ.get(env_var)
    if not candidate:
        return default
    path = Path(candidate).expanduser()
    if path.is_absolute():
        return path
    return (default.parent / path).resolve()


def build_pipeline_paths(script_file: Path) -> PipelinePaths:
    """Compute helper script locations, honoring environment overrides."""

    utils_root = _find_util_root(script_file)

    merge_default = utils_root / "merge_mrc_s_in_this_list_noise.py"
    if not merge_default.exists():
        merge_default = (
            utils_root
            / "prepare_data"
            / "mrc"
            / "1_prepare_2d"
            / "w_cistem"
            / "2_mrc_files_mrcs"
            / "3_mrc_files_to_mrcs"
            / "merge_mrc_s_in_this_list_noise.py"
        )

    phase_template_default = utils_root / "in_nm.defocus_template"
    if not phase_template_default.exists():
        phase_template_default = (
            utils_root
            / "prepare_data"
            / "mrc"
            / "1_prepare_2d"
            / "w_cistem"
            / "3_phaseflip"
            / "in_nm.defocus_template"
        )

    reconstruct_default = utils_root / "2_reconstruct_simulated.bash"
    if not reconstruct_default.exists():
        reconstruct_default = (
            utils_root
            / "prepare_data"
            / "mrc"
            / "2_reconstruct_from_2D"
            / "2_reconstruct_simulated.bash"
        )

    return PipelinePaths(
        scripts_root=utils_root,
        merge_script=_resolve_path("GRIPTOMO_MERGE_SCRIPT", merge_default),
        phase_template=_resolve_path("GRIPTOMO_PHASE_TEMPLATE", phase_template_default),
        reconstruct_script=_resolve_path(
            "GRIPTOMO_RECONSTRUCT_SCRIPT", reconstruct_default
        ),
    )


def main() -> None:
    """Entrypoint for single-folder reconstruction."""

    start_time = time.time()
    start_dir = Path.cwd()

    if len(sys.argv) < 5:
        print(
            "Specify <min_tilt> <max_tilt> <increment> <low_pass_level> (e.g. python at_single_folder.py -60 60 3 0.16)"
        )
        raise SystemExit(1)

    min_tilt = int(sys.argv[1])
    max_tilt = int(sys.argv[2])
    increment = float(sys.argv[3])
    low_pass_level = float(sys.argv[4])

    python_bin = os.environ.get("GRIPTOMO_PYTHON_BIN", sys.executable)
    ctfphaseflip_bin = (
        os.environ.get("IMOD_CTFPHASEFLIP_BIN")
        or os.environ.get("GRIPTOMO_CTFPHASEFLIP_BIN")
        or "ctfphaseflip"
    )

    paths = build_pipeline_paths(Path(__file__).resolve())

    angles = generate_angles(min_tilt, max_tilt, increment)

    if min_tilt < 0:
        rec_dir_name = f"m{abs(min_tilt)}_{max_tilt}_{increment}"
    else:
        rec_dir_name = f"{min_tilt}_{max_tilt}_{increment}"
    rec_dir = start_dir / rec_dir_name

    if rec_dir.exists():
        print(f"{rec_dir_name} exists already. Exiting without changes.")
        raise SystemExit(1)
    rec_dir.mkdir()

    restore_previous_run(start_dir)

    mrc_to_mrcs(start_dir, angles, paths, python_bin)
    phase_flip(start_dir, angles, paths, ctfphaseflip_bin)
    reconstruct_from_2d(start_dir, paths)

    move_outputs(start_dir, rec_dir)
    apply_filters(rec_dir, low_pass_level)
    archive_remaining(start_dir)

    duration = time.time() - start_time
    print(f"reconstruction took {duration:.2f}s")


if __name__ == "__main__":
    main()
