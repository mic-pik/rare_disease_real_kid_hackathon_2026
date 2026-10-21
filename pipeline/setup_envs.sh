#!/usr/bin/env bash
# Install recorded direct dependencies. Requires uv and a compatible CUDA toolkit.
set -euo pipefail
TASK_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$TASK_ROOT"
uv venv -p 3.13 track1/.venv
uv pip install --python track1/.venv/bin/python -r envs/esm.txt
uv venv -p 3.12 track1/.venv-evo2
uv pip install --python track1/.venv-evo2/bin/python 'torch==2.7.1+cu128' --index-url https://download.pytorch.org/whl/cu128
uv pip install --python track1/.venv-evo2/bin/python setuptools wheel packaging ninja psutil numpy
# flash-attn needs torch available at build time. Set CUDA_HOME to your toolkit if necessary.
uv pip install --python track1/.venv-evo2/bin/python 'flash-attn==2.8.0.post2' --no-build-isolation
uv pip install --python track1/.venv-evo2/bin/python -r envs/evo.txt --extra-index-url https://download.pytorch.org/whl/cu128
uv venv -p 3.11 gene_mol/.venv
uv pip install --python gene_mol/.venv/bin/python -r envs/chem.txt
uv venv -p 3.12 mva_llm/.venv-vllm
uv pip install --python mva_llm/.venv-vllm/bin/python -r envs/reasoning.txt
for env in track1/.venv track1/.venv-evo2 gene_mol/.venv mva_llm/.venv-vllm; do
    uv pip freeze --python "$env/bin/python" > "$env/resolved-requirements.txt"
done
