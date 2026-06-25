import os
import sys
from pathlib import Path

HELPER_BIN = os.environ.get("CISTEM_HELPER_BIN", "python3")
ANGLES_PATH = Path(__file__).with_name("angles.txt")

try:
    ANGLES = [
        float(line.strip())
        for line in ANGLES_PATH.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
except FileNotFoundError as exc:
    raise SystemExit(f"angles.txt is required beside {__file__}") from exc

for doseval in [replace_this]:
    # doseval ('x') is electrons per Ang^2 per frame; keep <= 20 and typically use 1-2 to stay near 100-120 total dose.
    print(f"Working on dose = {doseval} now")

    for frameval in [replace_this]:
        # frameval ('e') counts frames per tilt; doseval * frameval yields per-tilt dose and scales series dose.
        print(f"Working on frame = {frameval} now")

        for thickval in [replace_this]:
            # thickval ('t') is ice thickness in Angstroms; lamella runs hover near 3000 and values < 1 fail cisTEM.
            print(f"Working on thickness = {thickval} now")

            for phaseval in [replace_this]:
                # phaseval models the phase plate offset; 0.5 reflects typical phase-plate acquisition while 0 enables no-plate controls.
                print(f"Working on extra phase = {phaseval} now")

                # Example: doseval 20 with frameval 1 across 361 tilts (~0.5 deg) yields ~7240 electrons per Ang^2, far above experimental norms.

                new_dir = f"simulate_{doseval}x{frameval}e_{thickval}t_{phaseval}p"
                if not os.path.isdir(new_dir):
                    os.mkdir(new_dir)
                os.chdir(new_dir)

                for angleval in ANGLES:
                    # Legacy coarse 3 degree angles retained for reference: [-60.0,-57.0,...,60.0]
                    print(f"Working on angle = {angleval} now")

                    for defocusval in [replace_this]:
                        # defocusval should be >= 15000 Angstroms (~1.5 micron) for typical tilt series.
                        print(f"Working on defocus = {defocusval} now")

                        print(f"Replace defocusval of star file with {defocusval} now")

                        if not os.path.exists("../replace_defocus_in_star.py"):
                            print("replace_defocus_in_star.py is not found.")
                            print(
                                "Copy replace_defocus_in_star.py into this directory and retry."
                            )
                            sys.exit(1)
                        cmd = f"{HELPER_BIN} ../replace_defocus_in_star.py {defocusval}"
                        print(cmd)
                        os.system(cmd)

                        cmd = f"source ../simulate-tilt-noise.sh {doseval} {frameval} {thickval} {phaseval} {angleval} {defocusval}"
                        print(cmd)
                        os.system(cmd)

                os.chdir("..")
