# Parsl HPC Quickstart

This guide outlines the bare essentials for running GRIP-Tomo 2.0 workflows with Parsl on a generic HPC system. Adapt the examples to match your site's scheduler, module layout, and security policies.

## 1. Prerequisites
- Install GRIP-Tomo with Parsl extras in an isolated environment:
  ```bash
  uv venv
  source .venv/bin/activate
  uv pip install .[parsl]
  ```
- Confirm remote nodes can access the same Python environment (shared filesystem or site module).
- Collect scheduler details: queue/partition name, walltime limits, account/allocation code, node and core counts.

## 2. Configure Parsl
Parsl uses a `Config` object to describe how workers launch. Copy `griptomo/parsl/parsl_config_ci.py` as a template and update the provider section for your cluster. A minimal Slurm example:

```python
from parsl.addresses import address_by_hostname
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.providers import SlurmProvider

config = Config(
    executors=[
        HighThroughputExecutor(
            label="griptomo-htex",
            address=address_by_hostname(),
            max_workers=1,
            provider=SlurmProvider(
                partition="gpu",
                nodes_per_block=1,
                cores_per_node=32,
                walltime="02:00:00",
                init_blocks=1,
                max_blocks=4,
                scheduler_options="#SBATCH --account=my_allocation",
                worker_init="source /path/to/.venv/bin/activate",
            ),
        )
    ],
)
```

Key fields to adjust:
- `partition`, `scheduler_options`: match your queue and allocation.
- `nodes_per_block`, `cores_per_node`, `max_blocks`: tune for available resources.
- `worker_init`: load modules, activate the virtual environment, or set environment variables.

## 3. Verify Connectivity
Before launching the full pipeline, run a quick Parsl sanity check on the login node:
```bash
uv run python - <<'PY'
from parsl import load, python_app
from griptomo.parsl import parsl_config_ci

load(parsl_config_ci.config)

@python_app
def whoami():
    import socket
    return socket.gethostname()

print(whoami().result())
PY
```
If the job completes and prints a worker hostname, the configuration is valid. Troubleshoot common issues (authentication, walltime, module loads) before proceeding.

## 4. Launch GRIP-Tomo Tasks
1. Stage input data on a shared filesystem accessible to compute nodes.
2. Ensure required third-party binaries (IMOD, EMAN2, RELION, cisTEM) are available in `worker_init`.
3. Submit your Parsl-enabled script (for example `griptomo/parsl/apps.py`) with:
   ```bash
   uv run python griptomo/parsl/apps.py --config griptomo.parsl.parsl_config_ci
   ```
   Adjust CLI parameters to match your datasets and output locations.

## 5. Operational Tips
- **Logging**: Parsl logs to `runinfo/`; inspect the newest directory when debugging.
- **Scaling**: Start with a single block (`init_blocks=1`, `max_blocks=1`) and increase after validation.
- **Retries**: Enable Parsl retries (`retries` parameter on apps) for long-running workloads.
- **Profiling**: Set `monitoring=True` in the config to collect execution metrics if your site allows it.

## 6. Security & Policy Checklist
- Verify that automated job submission complies with site policies.
- Use site-provided secure storage for credentials or SSH keys.
- Clean up temporary data and log directories when runs complete.

With these steps you should be able to adapt the Parsl configuration to most schedulers. For advanced patterns (MPI workloads, heterogeneous resources), refer to the Parsl docs: https://parsl.readthedocs.io.
