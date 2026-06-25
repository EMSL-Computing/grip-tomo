# AI Agents & LLM-Assisted Development Guide

This document provides guidance for AI coding agents (GitHub Copilot, Cursor, Windsurf, Cline, Claude, etc.) and human developers using LLMs to contribute to GRIP-Tomo 2.0.

## Quick Start for AI Agents

**Primary Instructions**: See [`.github/copilot-instructions.md`](.github/copilot-instructions.md) for comprehensive architectural guidance, coding patterns, and conventions.

**Key Resources**:
- **Core Pipeline**: `griptomo/core/{density2graph,pdb2graph,graph2class}.py`
- **Test Examples**: `tests/test_core.py` demonstrates all major APIs
- **Tutorial**: `docs/tutorial_notebook/tutorial.ipynb` shows complete workflows
- **Utils Scripts**: `utils/*.py` for simulation orchestration patterns

## AI-Friendly Patterns

### 1. Dual Backend Support (CPU/GPU)
Always check availability flags before dispatching:
```python
from griptomo.core import density2graph as d2g

if d2g.IS_GPU_AVAILABLE:
    # Use cuML/cuGraph implementations
    result = gpu_function(data)
else:
    # Fallback to CPU (scikit-learn/NetworkX)
    result = cpu_function(data)
```

### 2. Type Flexibility
Functions accept multiple input types via `_ensure_array()`:
```python
# All valid inputs:
mrc_file = mrcfile.mmap("data.mrc")
numpy_array = np.array([...])
result1 = d2g.normalize_mrc_data(mrc_file)      # MrcFile object
result2 = d2g.normalize_mrc_data(numpy_array)   # numpy array
result3 = d2g.normalize_mrc_data(mrc_file.data) # data attribute
```

### 3. Parsl Self-Contained Apps
All imports must be inside the function body:
```python
from parsl import python_app

@python_app
def process_data(input_path):
    # ✅ Correct: imports inside function
    import numpy as np
    from griptomo.core import density2graph as d2g
    
    # Process data...
    return result

# ❌ Wrong: imports at module level won't work on remote workers
```

### 4. Environment Variable Patterns
Always provide fallback defaults:
```python
import os
from pathlib import Path

# Standard pattern for tool paths
TOOL_BIN = os.getenv("GRIPTOMO_TOOL_BIN", "default_tool_name")

# Directory paths should use Path objects
OUTPUT_DIR = Path(os.getenv("GRIPTOMO_OUTPUT_DIR", "./output"))
```

### 5. Shell Command Construction
Use `shlex.join()` for safe command building:
```python
import shlex
import subprocess

cmd = [tool_bin, "--input", str(input_path), "--output", str(output_path)]
print(f"Running: {shlex.join(cmd)}")
result = subprocess.run(cmd, capture_output=True, text=True)
```

## Common Code Generation Scenarios

### Adding a New Feature Extraction Function

**Pattern to follow** (`griptomo/core/graph2class.py`):
```python
def extract_new_feature(G: nx.Graph) -> float:
    """
    Brief description of the feature.
    
    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    
    Returns
    -------
    float
        The computed feature value.
    
    Notes
    -----
    Describe algorithm, computational complexity, or references.
    """
    # Implementation
    return feature_value
```

**Add corresponding test** (`tests/test_core.py`):
```python
def test_new_feature(self):
    """Tests new feature extraction."""
    G = nx.karate_club_graph()  # Or use fixture data
    result = g2c.extract_new_feature(G)
    assert isinstance(result, float)
    assert result > 0  # Add meaningful assertions
```

### Adding a Utility Script

**Template structure**:
```python
#!/usr/bin/env python
"""Brief description of what this script does."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

def main(argv=None):
    """Main entry point with CLI argument parsing."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input file/directory")
    parser.add_argument("--output", type=Path, default=Path("output"))
    args = parser.parse_args(argv)
    
    # Respect environment variables
    tool_bin = os.getenv("GRIPTOMO_TOOL_BIN", "default_tool")
    
    # Implementation...
    print(f"Processing {args.input}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
```

**Add corresponding test** (`tests/test_new_script.py`):
```python
import pytest
from pathlib import Path

def test_new_script_smoke(tmp_path, monkeypatch):
    """Basic smoke test for new script."""
    # Import the module
    from utils import new_script
    
    # Create test inputs
    test_input = tmp_path / "input.txt"
    test_input.write_text("test data")
    
    # Override environment
    monkeypatch.setenv("GRIPTOMO_TOOL_BIN", "/mock/tool")
    
    # Run and assert
    result = new_script.main([str(test_input)])
    assert result == 0
```

## Testing Guidance

### Run Tests Before Committing
```bash
# Quick smoke test (core modules only)
pytest -q tests/test_core.py

# Full suite
pytest -v

# Specific test file
pytest tests/test_density2graph.py -v

# Skip GPU tests if not available
pytest -v -m "not gpu"
```

### Test Patterns to Follow
1. **Use fixtures**: `tmp_path` for temporary directories, `monkeypatch` for environment isolation
2. **Mock external binaries**: Use `monkeypatch.setattr()` to intercept subprocess calls
3. **Parametrize**: Test multiple backends/configurations in one test function
4. **Integration tests**: Mark long-running tests with `@pytest.mark.slow` or `@pytest.mark.integration`

## Documentation Standards

### Docstring Format (NumPy Style)
```python
def function_name(param1: type1, param2: type2) -> return_type:
    """
    One-line summary (imperative mood).
    
    Longer description if needed. Can span multiple paragraphs.
    
    Parameters
    ----------
    param1 : type1
        Description of param1.
    param2 : type2
        Description of param2.
    
    Returns
    -------
    return_type
        Description of return value.
    
    Raises
    ------
    ValueError
        When invalid input is provided.
    
    Notes
    -----
    Additional context, algorithms, references, or warnings.
    
    Examples
    --------
    >>> result = function_name(1, 2)
    >>> print(result)
    3
    """
```

### Regenerate API Docs After Changes
```bash
uv run pdoc griptomo --docformat numpy --output-directory docs/api
```

## Known Edge Cases & Gotchas

1. **MRC File Memory Mapping**: Always use `mode='r'` for read-only access. Never write to memory-mapped files.
   
2. **Graph Connectivity**: Many feature extraction functions require connected graphs. Use `nx.is_connected()` and extract largest component if needed:
   ```python
   if not nx.is_connected(G):
       G_cc = sorted(nx.connected_components(G), key=len, reverse=True)
       G = G.subgraph(G_cc[0])
   ```

3. **Parsl Execution Context**: Parsl apps run on remote workers without access to the parent process's environment. All dependencies must be explicitly imported inside the function.

4. **Platform-Specific GPU Packages**: GPU extras only work on Linux x86_64. The code gracefully degrades to CPU on other platforms, but tests should check `IS_GPU_AVAILABLE` flag.

5. **Path Handling**: Use `pathlib.Path` for cross-platform compatibility. Convert to `str` only when passing to external subprocess commands.

## Development Workflow

1. **Create feature branch**: `git checkout -b feature/new-feature`
2. **Write code + tests**: Follow patterns above
3. **Run tests**: `pytest -v`
4. **Update docs if needed**: Regenerate API docs, update quickstarts
5. **Commit with descriptive message**: Follow conventional commits style
6. **Push and create PR**: Target `griptomo_2_prerelease` branch

## Getting Help

- **Architecture questions**: See `.github/copilot-instructions.md`
- **API usage**: Check `docs/tutorial_notebook/tutorial.ipynb`
- **Workflow specifics**: See `docs/cistem_subtomogram_simulation_quickstart.md` or `docs/parsl_hpc_quickstart.md`
- **Test examples**: Browse `tests/` directory for patterns

## Agent Performance Tips

When generating code for GRIP-Tomo:
1. **Always check existing patterns first** - search for similar functions in the codebase
2. **Prefer composition over duplication** - reuse existing utility functions
3. **Test incrementally** - run pytest after each change to catch issues early
4. **Respect conventions** - follow the project's import style, docstring format, and naming patterns
5. **Document assumptions** - if you're uncertain about behavior, add a comment or TODO

---

**Last Updated**: November 18, 2024  
**For AI Agents**: This file is intended to be read by both LLMs and humans. When making code suggestions, reference specific patterns from this guide and `.github/copilot-instructions.md`.
