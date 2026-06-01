# GRIP-Tomo 2.0.0 Alpha Release Notes

This document describes the changes introduced in GRIP-Tomo 2.0.0-alpha compared to version 1.0.

**Release Date**: June 1 2026 (v2), November 18, 2024 (v1)

### Added
- **GPU Acceleration**: Optional RAPIDS (cuDF, cuML, cuGraph) support for clustering and feature extraction on NVIDIA GPUs (Linux x86_64 only)
- **Parsl HPC Workflows**: Distributed computing support for HPC clusters with Slurm/PBS via `griptomo.parsl` module
- **Comprehensive Test Suite**: 19 test modules with 70+ test functions covering core modules, utils, and integration scenarios
- **API Documentation**: Auto-generated pdoc documentation in `docs/api/` with NumPy-style docstrings
- **Tutorial Notebook**: Interactive Jupyter notebook (`docs/tutorial_notebook/tutorial.ipynb`) with complete workflow examples
- **Quickstart Guides**: Three detailed guides for subtomogram simulation (cisTEM), Parsl HPC setup, and general usage
- **Environment Configuration**: `.env.example` template for managing external tool paths and runtime settings
- **igraph Backend**: CPU-optimized feature extraction via `ig_extract_features.py` alongside NetworkX implementation
- **Modern Packaging**: Migration to `pyproject.toml` with uv-based environment management and pinned dependencies
- **Utils Pipeline**: Complete set of orchestration scripts for cisTEM simulation workflows (`utils/`)
  - `prepare_multiple_conditions.py`: Batch input generation for parameter sweeps
  - `at_single_folder.py`: Post-cisTEM reconstruction pipeline (merge, CTF, filtering, inversion)
  - `use_parallelized_cores.py`: Parallel execution driver for multi-folder processing
- **AI Coding Support**: Agent instruction files (`.github/copilot-instructions.md`, `AGENTS.md`) for LLM-assisted development

### Changed
- **Breaking**: Refactored module structure from flat layout to `griptomo.core` and `griptomo.parsl` packages
- **Breaking**: Deprecated `normalize_and_threshold_data()` in favor of modular `normalize_mrc_data()` + `threshold_mrc_data()` + `generate_point_cloud_from_mrc_data()` pipeline
- **Breaking**: Minimum Python version raised to 3.10
- **Dependency Updates**: numpy 1.26.4, networkx 3.3, pandas 2.2.2, scikit-learn 1.5.1, scipy 1.14.0, igraph 0.11.6
- **Documentation**: Switched from Sphinx to pdoc for API documentation generation
- **Type Annotations**: Added `from __future__ import annotations` for forward compatibility and better type hinting
- **File Handling**: Improved MRC file memory-mapping with explicit read-only mode for large files

### Fixed
- Repository name typo in README citations
- MRC data normalization edge cases with explicit bounds checking
- Graph connectivity validation in feature extraction functions

### Removed
- **Breaking**: Sphinx-based documentation build system and `docs/_build/` artifacts
- Python 2.x compatibility code and legacy string formatting
- ReadTheDocs configuration (`.readthedocs.yaml`)
- Old GitHub Actions workflow (`python_build.yml`) - replaced by local test execution

### Known Issues
- `leave_2D_density()` function lacks complete documentation (TODO at line 813 in `density2graph.py`)
- `test_identify_threshold_ratio()` unit test needs review/validation (TODO at line 715 in `test_core.py`)

### Migration Notes
**From v1.0 to v2.0:**
1. Update imports: `from griptomo.core import density2graph` instead of `import density2graph`
2. Replace deprecated function:
   ```python
   # Old (v1.0)
   xyz_data = d2g.normalize_and_threshold_data(mrc, threshold)
   
   # New (v2.0)
   D_norm = d2g.normalize_mrc_data(mrc.data)
   D_thresh = d2g.threshold_mrc_data(D_norm, threshold)
   xyz_data = d2g.generate_point_cloud_from_mrc_data(D_thresh, threshold)
   ```
3. Install with extras: `uv pip install -e ".[dev]"` for development, add `[parsl]` for HPC, `[gpu-cu12]` for GPU
4. Review `.env.example` and configure external tool paths for cisTEM workflows
