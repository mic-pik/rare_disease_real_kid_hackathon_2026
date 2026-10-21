#!/usr/bin/env bash
# Submit from the repository root; supply account/partition/resources to sbatch.
# Example: sbatch --gres=gpu:4 --cpus-per-task=32 --mem=300G --time=2-00:00:00 pipeline/run.sbatch --config config/local.json --stage all
set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Submit this script with sbatch from the repository root}"
cd "$SLURM_SUBMIT_DIR"
export VLLM_WORKER_MULTIPROC_METHOD=spawn
python3 pipeline/run.py "$@"
