try:
    import pretty_errors  # type: ignore[unused-import]
except Exception:
    print("pretty-errors not installed; continuing without enhanced tracebacks.")

import os
import sys

# run @each_input_folders

starting_dir = os.getcwd()


def submit_at_many_folders():
    sbatch_bin = os.environ.get("SLURM_SBATCH", "sbatch") or "sbatch"
    for folder_candidate in os.listdir(starting_dir):
        each_pdb_input_folder = os.path.join(starting_dir, folder_candidate)
        print(f"each_pdb_input_folder:{each_pdb_input_folder}")
        if os.path.isdir(each_pdb_input_folder):
            for each_condition_folder in os.listdir(each_pdb_input_folder):
                print(f"each_condition_folder:{each_condition_folder}")

                if each_condition_folder == "stars":
                    print("Run submit_at_many_folders.py at each_input_folders")
                    sys.exit(1)

                if "template_input" not in each_condition_folder:
                    sbatch_having_folder = os.path.join(
                        each_pdb_input_folder, each_condition_folder
                    )
                    print(f"sbatch_having_folder:{sbatch_having_folder}")

                    # check whether sbatch_having_folder is a folder/directory
                    if not os.path.isdir(sbatch_having_folder):
                        print(f"{sbatch_having_folder} is not a folder/directory")
                        sys.exit(1)

                    original_cwd = os.getcwd()
                    try:
                        os.chdir(sbatch_having_folder)
                        cmd = f"{sbatch_bin} run.sbatch"
                        # print(cmd)
                        result = os.system(cmd)
                        if result != 0:
                            print(
                                f"Submission command failed with exit code {result} for {sbatch_having_folder}"
                            )
                    finally:
                        os.chdir(original_cwd)


###### end of def submit_at_many_folders():

if __name__ == "__main__":
    if not hasattr(sys, "version_info") or sys.version_info < (3, 7):
        raise ValueError("Script requires Python 3.7 or higher!")

    submit_at_many_folders()
