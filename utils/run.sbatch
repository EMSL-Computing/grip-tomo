#!/usr/bin/env bash

#SBATCH -A PROJECT_ACCOUNT            # TODO: replace with your allocation/account ID
#SBATCH --time=02-00:00:00            # TODO: adjust walltime (DD-HH:MM:SS)
#SBATCH --nodes=1                     # TODO: set node count as required
#SBATCH --job-name=cisTEM             # Optional: rename the job
#SBATCH --output=cisTEM.%j.out        # STDOUT log path
#SBATCH --error=cisTEM.%j.err         # STDERR log path
#SBATCH --partition=compute           # TODO: pick the appropriate queue/partition
##SBATCH --gpus=1                     # Uncomment if GPUs are required

# --- Environment setup (customize for your HPC site) ---
# module purge
# module load cisTEM/2.0.0-alpha
# source /path/to/conda_or_virtualenv/bin/activate
# export CISTEM_BIN=/cluster/apps/cisTEM/2.0.0-alpha

python runscript.py
