import os
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _write_executable(path: Path, content: str) -> None:
    path.write_text(content)
    path.chmod(0o755)


def test_reconstruct_simulated_bash_smoke(tmp_path):
    """Smoke test for the 2_reconstruct_simulated.bash orchestration with stubbed IMOD binaries."""
    script_path = REPO_ROOT / "utils" / "2_reconstruct_simulated.bash"

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()

    tilt_log = tmp_path / "tilt.log"

    _write_executable(
        bin_dir / "header",
        "#!/bin/bash\necho 'Number of columns, rows, sections : 64 64 1'\n",
    )

    _write_executable(
        bin_dir / "tilt",
        "#!/bin/bash\n"
        f"printf '%s ' \"$@\" >> '{tilt_log}'\n"
        f"printf '\n' >> '{tilt_log}'\n"
        "output=''\n"
        "while [[ $# -gt 0 ]]; do\n"
        "  case $1 in\n"
        "    -output) shift; output=$1 ;;\n"
        "  esac\n"
        "  shift\n"
        "done\n"
        'if [[ -n "$output" ]]; then\n'
        '  touch "$output"\n'
        "fi\n",
    )

    _write_executable(
        bin_dir / "trimvol",
        "#!/bin/bash\nprintf '%s\n' \"$@\" >> '${TMPDIR:-/tmp}/trimvol.log'\n",
    )

    sample = tmp_path / "sample.mrcs"
    sample.write_bytes(b"")

    (tmp_path / "angles.rawtlt").write_text("0\n")
    (tmp_path / "placeholder.rec~").write_bytes(b"")

    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["IMOD_TILT"] = str(bin_dir / "tilt")
    env["IMOD_TRIMVOL"] = str(bin_dir / "trimvol")

    result = subprocess.run(
        ["/bin/bash", str(script_path)],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "sample.mrcs.rec").exists()
    assert tilt_log.exists()
    assert "-input" in tilt_log.read_text()
