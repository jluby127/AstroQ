#!/usr/bin/env bash
# One-time HIRES ops setup. Run from repo root:
#   bash ops/HIRES/setup_env.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CONDA_ENV="${CONDA_ENV:-astroq}"

echo "Repo:       $REPO_ROOT"
echo "Conda env:  $CONDA_ENV"

conda run -n "$CONDA_ENV" pip install -e "$REPO_ROOT"
# HIRES-only packages (jump2-client for past.csv), from a private index.
conda run -n "$CONDA_ENV" pip install -U -r "$REPO_ROOT/ops/HIRES/requirements.txt"
conda env config vars set -n "$CONDA_ENV" CC_OUTPUT_PATH="$REPO_ROOT"

echo ""
echo "Done. Reactivate and verify:"
echo "  conda deactivate && conda activate $CONDA_ENV"
echo "  python -c \"import astroq; print(astroq.__file__)\""
echo "  python -c \"import jump2_client\""
echo "  echo \$CC_OUTPUT_PATH"
echo ""
echo "jump2 credentials (needed by prep to build past.csv): create an API key"
echo "under Profile -> API Keys on https://cannon.caltech.edu, then either"
echo "  export JUMP_API_KEY=<key>   (e.g. in ops/HIRES/.env)"
echo "or run  jump2 scripts login"
