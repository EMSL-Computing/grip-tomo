# Subtomogram Simulation Workflow Quickstart

This guide walks through preparing, running, and post-processing cisTEM simulations for multiple PDB structures and experimental conditions using the bundled `utils/` helpers. The end result is a set of reconstructed 3D subtomograms that have been filtered and density-inverted.

## 1. Prerequisites
- **Software**: Python 3.10+, cisTEM, IMOD (`tilt`, `trimvol`, `clip`, `ctfphaseflip`), EMAN2 (`e2proc2d.py`), RELION (`relion_image_handler`).
- **Environment**: clone the repository, create a uv-managed environment, and install GRIP-Tomo with any extras you need:
  ```bash
  uv venv
  source .venv/bin/activate
  uv pip install .           # add [parsl], [gpu-cu11], or [gpu-cu12] as needed
  # optional for active development
  uv pip install -e .
  ```
- **Data**: a folder of input `.pdb` files and a conditions file patterned after `utils/try_these_conditions.txt`.

## 2. Environment Configuration
The workflow relies on environment variables to locate external tools. A template is provided at `.env.example` in the repository root.

1. Copy `.env.example` to `.env` and edit the required paths:
   ```bash
   cp .env.example .env
   # Edit .env with your system's paths for:
   # - GRIPTOMO_OUTPUT_DIR (base output directory)
   # - External tool binaries (IMOD_CTFPHASEFLIP_BIN, IMOD_TILT, EMAN2_E2PROC2D, etc.)
   # - cisTEM installation paths
   ```
2. Source the file before running any helpers:
   ```bash
   source .env
   ```
3. Optional: set `FORCE=1` to allow regenerating inputs if `each_input_folders/` already exists.

## 3. Prepare Simulation Inputs
Run the preparation script from the repository root:
```bash
python utils/prepare_multiple_conditions.py \
  path/to/try_these_conditions.txt \
  path/to/folder_with_pdbs
```
This creates `each_input_folders/<pdb>/<condition>/template_input/` populated with the necessary run scripts, stars, and defocus templates.

## 4. Launch cisTEM Simulations
From within `each_input_folders/`, either launch Slurm jobs via the provided template or run cisTEM locally.

- **Slurm example** (uses `SLURM_SBATCH` if exported):
  ```bash
  cd each_input_folders
  python ../utils/submit_at_many_folders.py
  ```
- **Manual run**: enter a condition directory and execute `runscript.py` / `simulate-tilt-noise.sh` as needed.

Wait for all cisTEM runs to finish; each condition directory should contain `cisTEM.*` outputs and `simulate_*` folders.

## 5. Collect Outputs Into `output/`
After simulations finish, normalize the directory layout:
```bash
cd each_input_folders
python ../utils/run_at_many_input_folders_that_has_both_many_pdb_and_many_conditions.py
```
This moves cisTEM artefacts into `output/` under every condition so the reconstruction step can find them.

## 6. Reconstruct, Filter, and Invert
Every `output/simulate_*` directory must run the single-folder pipeline. There are two common approaches:

### Option A – Parallel driver
Execute once from `each_input_folders/` after activating your environment:
```bash
cd each_input_folders
python - <<'PY'
from utils import use_parallelized_cores
use_parallelized_cores.starting_dir = '.'
use_parallelized_cores.modified_run_at_many_folders(
    run_this_py='../../utils/at_single_folder.py',
    min_tilt=-60,
    max_tilt=60,
    increment=3,
    low_pass_level=0.16,
)
PY
```
Adjust the tilt range and low-pass level for your experiment. The helper honors `GRIPTOMO_PYTHON_BIN` and `GRIPTOMO_MAX_WORKERS` when set.

### Option B – Manual loop
If you prefer shell loops, run `at_single_folder.py` inside each simulation directory:
```bash
find each_input_folders -type d -name "simulate_*" -print0 | \
while IFS= read -r -d '' sim_dir; do
  (cd "$sim_dir" && python ../../../../utils/at_single_folder.py -60 60 3 0.16)
done
```

During execution the script will:
1. Merge noise-projected `.mrc` files into an `.mrcs` tilt stack (`merge_mrc_s_in_this_list_noise.py`).
2. Apply CTF phase flipping via `ctfphaseflip` using `in_nm.defocus_template`.
3. Reconstruct 3D volumes with IMOD (`2_reconstruct_simulated.bash`).
4. Low-pass filter and density-invert the resulting `.rec` files (`filter_mrc.py`, `invert_mrc_density.py`).

Each simulation produces a reconstruction folder named after the tilt range (for example, `m60_60_3`) containing the filtered and inverted tomograms.

## 7. Where to Look Next
- Final volumes: `each_input_folders/<pdb>/<condition>/output/simulate_*/m60_60_3/*.rec` (or the name matching your tilt parameters).
- Logs and intermediate CSVs remain in the same reconstruction folder for traceability.
- Rerun any step by deleting the target reconstruction directory and setting `FORCE=1` (for input regeneration) if needed.

With these steps you can batch-process multiple PDB structures across parameter grids and obtain filtered, inverted tomograms ready for downstream analysis.
