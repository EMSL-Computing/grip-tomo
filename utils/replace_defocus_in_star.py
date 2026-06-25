from __future__ import annotations

import glob
import os
import shutil
import sys
from pathlib import Path


STAR_ROOT = os.environ.get("STAR_DIR", "../stars")


def replace_defocus_in_star(defocus_in_angstrom: int, pattern: str) -> None:
    for input_star in glob.glob(pattern):
        print(f"processing {input_star}")
        input_path = Path(input_star)
        output_path = input_path.with_name(f"{input_path.stem}_new.star")

        with input_path.open("r", encoding="utf-8") as ref_star_file:
            ref_star_file_lines = ref_star_file.readlines()

        with output_path.open("w", encoding="utf-8") as output_star_file:
            pos_found = False
            for line in ref_star_file_lines:
                if not pos_found:
                    output_star_file.write(line)
                    if line[:8].rstrip() == "#    POS":
                        pos_found = True
                else:
                    split_line = line.split()
                    new_fields = []
                    for index, value in enumerate(split_line):
                        if index in (6, 7):
                            new_fields.append(str(defocus_in_angstrom))
                        else:
                            new_fields.append(value)
                    output_star_file.write("   " + " ".join(new_fields) + "\n")

        shutil.move(str(output_path), input_star)


if __name__ == "__main__":
    args = sys.argv[:]

    if len(args) < 2:
        print("usage: python replace_defocus_in_star.py <defocus_angstrom>")
        sys.exit(1)

    defocus_in_angstrom = int(args[1])
    glob_pattern = os.path.join(STAR_ROOT, "tilt_*.star")
    replace_defocus_in_star(defocus_in_angstrom, glob_pattern)
