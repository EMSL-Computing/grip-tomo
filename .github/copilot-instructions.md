# GRIP-Tomo 2.0 AI Coding Agent Instructions

## Project Overview
GRIP-Tomo converts cryo-electron tomography (cryo-ET) subtomograms into graphs, extracts features, and trains ML classifiers to identify proteins in 3D cellular volumes. This is research software under active development with potential for breaking changes.

**Core workflow**: Density maps (`.mrc` files) → Point clouds → Clustered centroids → Graphs → Feature extraction → ML classification

## Architecture & Key Modules

### Core Pipeline (`griptomo/core/`)
- **`density2graph.py`**: Converts 3D density maps to graphs via clustering (DBSCAN/HDBSCAN). Supports both CPU (scikit-learn) and GPU (RAPIDS cuML/cuGraph) backends. Check `IS_GPU_AVAILABLE` flag before using GPU features.
- **`pdb2graph.py`**: Generates reference graphs from PDB/PDBx protein structures. Uses BioPython for parsing and assigns edges based on distance cutoffs (default 8Å between alpha carbons).
- **`graph2class.py`**: Extracts graph-theoretic features (centrality, clustering, path lengths) for ML classification. Supports both NetworkX and igraph backends.
- **`ig_extract_features.py`** / **`gpu_extract_features.py`**: Backend-specific feature extraction optimized for CPU (igraph) or GPU (cuGraph).

**Key pattern**: Functions accept both file paths and in-memory objects (mrcfile.MrcFile, numpy/cupy arrays, networkx/igraph graphs). Use `_ensure_array()` for type normalization.

### Parsl HPC Workflows (`griptomo/parsl/`)
- **`apps.py`**: Defines `@python_app` decorated functions for distributed execution on HPC clusters (Slurm, PBS, etc.).
- **`parsl_config_ci.py`**: Template configuration for cluster-specific settings (partition, walltime, worker initialization).
- Parsl apps are self-contained: all imports must be inside the function body to work on remote workers.

### Utilities (`utils/`)
Scripts orchestrate multi-step simulations using external tools (cisTEM, IMOD, EMAN2, RELION):
- **`prepare_multiple_conditions.py`**: Batch-generates input directories for parameter sweeps across PDB structures and experimental conditions.
- **`at_single_folder.py`**: Post-cisTEM pipeline (merge tilt stacks, CTF correction, reconstruction, filtering, density inversion) for one simulation folder.
- **`use_parallelized_cores.py`**: Parallel driver that launches `at_single_folder.py` across multiple simulation outputs using `multiprocessing`.
- Scripts honor environment variables (`GRIPTOMO_*` prefixes) for binary paths, output directories, and parallelization settings.

**Convention**: Shell commands are constructed with `shlex.join()` and executed via `subprocess.run()`. Always capture and log both stdout and stderr.

## Environment & Dependencies

### Setup Commands
```bash
uv venv                                    # Create isolated environment
source .venv/bin/activate                   # Activate environment
uv pip install -e ".[dev]"                 # Dev install with testing tools
uv pip install -e ".[dev,parsl]"           # Add Parsl for HPC workflows
uv pip install -e ".[dev,gpu-cu12]"        # GPU support (Linux x86_64 + CUDA 12 only)
```

### Critical Dependencies
- **CPU baseline**: numpy 1.26.4, networkx 3.3, igraph 0.11.6, hdbscan 0.8.39, scikit-learn 1.5.1
- **GPU stack** (optional): cudf/cuml/cugraph 25.2.* from `https://pypi.nvidia.com` (Linux x86_64 only)
- **External tools**: IMOD (tilt, trimvol, clip, ctfphaseflip), EMAN2 (e2proc2d.py), RELION (relion_image_handler), cisTEM

GPU packages require NVIDIA GPUs and are platform-restricted. The codebase gracefully degrades to CPU when GPU is unavailable.

## Testing & Validation

### Run Tests
```bash
pytest -q tests/test_core.py              # Smoke tests for core modules
uv run pytest                              # Full test suite
pytest tests/test_parsl_apps.py           # Parsl workflow validation
```

### Test Structure
- **Fixtures**: `tests/fixtures/example_data/` contains minimal PDB/MRC files for unit tests.
- **Patterns**: Tests use `tmp_path` fixtures and `monkeypatch` for environment isolation. External binaries are mocked unless integration testing.
- **Naming**: `test_<module>_<behavior>` (e.g., `test_density2graph_converts_mrc_to_graph`).

Integration tests (`test_end_to_end_full_pipeline.py`) validate the complete simulation-to-reconstruction workflow but require external tools to be installed.

## Development Workflows

### Code Regeneration
```bash
uv run pdoc griptomo --docformat numpy --output-directory docs/api   # Regenerate API docs
```

### Adding New Features
1. **Core modules**: Add docstrings in NumPy format. Include `Parameters`, `Returns`, and `Notes` sections.
2. **Parsl apps**: Ensure all imports are inside the function body. Add `Returns` docstring describing the tuple structure.
3. **Utils scripts**: Accept CLI arguments via `argparse`. Respect `GRIPTOMO_*` environment variables for paths.
4. **Tests**: Add corresponding unit tests. Use parametrize for testing multiple backends (CPU/GPU, NetworkX/igraph).

### GPU-Specific Code
Always wrap GPU imports in try-except blocks and set fallback flags:
```python
try:
    import cupy as cp
    IS_GPU_AVAILABLE = cp.cuda.is_available()
except Exception:
    IS_GPU_AVAILABLE = False
```
Check `IS_GPU_AVAILABLE` before dispatching to GPU implementations. CPU code paths must work standalone.

## Key Conventions

### Data Normalization
- Density data is normalized to [-1, 1] via `normalize_mrc_data()` before thresholding.
- Threshold values are typically in range [0.0, 1.0] after normalization.
- **Inversion convention**: White voxels = high density. Some data sources require explicit inversion (`invert_mrc_density.py`).

### Graph Construction
- PDB graphs use alpha-carbon (CA) atoms by default (`CA_only=1`).
- Distance cutoffs (`d_cut`) typically 5-10 Angstroms for protein structure, variable for density-derived graphs.
- Graphs are stored as `.gexf` (NetworkX) or converted to igraph/cuGraph for feature extraction.

### File Organization
Simulation outputs follow structured hierarchy:
```
each_input_folders/<pdb>/<condition>/
  template_input/          # Generated by prepare_multiple_conditions.py
  output/simulate_*/       # cisTEM outputs, one per noise/tilt series
    m<min>_<max>_<inc>/   # Reconstruction directory (e.g., m60_60_3)
      *.rec               # Final reconstructed volumes
```

### Environment Variables
Scripts check these variables (fallback to defaults if unset):
- `GRIPTOMO_OUTPUT_DIR`: Base directory for simulation outputs
- `GRIPTOMO_PYTHON_BIN`: Python interpreter for subprocess calls (default: `sys.executable`)
- `GRIPTOMO_MAX_WORKERS`: Parallelization limit for `use_parallelized_cores.py`
- `IMOD_*`, `EMAN2_*`, `RELION_*`: Paths to external tool binaries

## Common Pitfalls

1. **MRC file modes**: Use `mrcfile.mmap(mode='r')` for read-only memory mapping on large files. Never write back to mapped files.
2. **Parsl retries**: Long-running workflows should enable `retries=N` on python_apps to handle transient cluster failures.
3. **NetworkX graph connectivity**: Many feature extraction functions require connected graphs. Use `nx.is_connected()` and extract largest component if needed (see `benchmark_igraph_graph_features()`).
4. **Path handling**: Use `pathlib.Path` for cross-platform compatibility. Convert to strings only when passing to external subprocess commands.
5. **Deprecated functions**: `normalize_and_threshold_data()` is deprecated; use `normalize_mrc_data()` + `threshold_mrc_data()` + `generate_point_cloud_from_mrc_data()` sequence.

## Documentation References
- Tutorial notebook: `docs/tutorial_notebook/tutorial.ipynb`
- Subtomogram simulation workflow: `docs/cistem_subtomogram_simulation_quickstart.md`
- Parsl HPC setup: `docs/parsl_hpc_quickstart.md`
- API docs (pdoc): `docs/api/index.html`

## Citation
When implementing features, preserve attribution comments (e.g., `# August George, 2022, PNNL`) at file headers. The project is published research software—see README for citation details.
