#!/bin/bash
# Build the conda environment TensorSignatures needs.
#
# WHY THIS IS ITS OWN ENVIRONMENT
# -------------------------------
# TensorSignatures 0.5.0 is pinned to a 2019 numerical stack - tensorflow <= 1.15,
# numpy < 1.17, scipy 1.3.2, scikit-learn 0.21.3, pandas 0.25.3 - and TensorFlow
# 1.15 was never built for Python >= 3.8: its wheels stop at cp37. It therefore
# CANNOT share an interpreter with a modern Python and must live in its own
# Python 3.7 environment.
#
# No root access is needed. Everything lands under $HOME and is removable with
#   rm -rf ~/miniconda3
# The system Python and R installations are untouched.
#
# Usage:  setup_tensorsig_env.sh
set -euo pipefail

CONDA_ROOT="${CONDA_ROOT:-$HOME/miniconda3}"
ENV_NAME="${ENV_NAME:-tensorsig}"
ENV_PY="$CONDA_ROOT/envs/$ENV_NAME/bin/python"

if [ -x "$ENV_PY" ] && "$ENV_PY" -c "import tensorsignatures" 2>/dev/null; then
  echo "environment already usable: $ENV_PY"
  "$ENV_PY" -c "import tensorsignatures as ts, tensorflow as tf, numpy as np; \
print('tensorsignatures', ts.__version__, '| tensorflow', tf.__version__, \
'| numpy', np.__version__)"
  exit 0
fi

# ---- miniconda pinned to the last py37 installer -----------------------------
if [ ! -x "$CONDA_ROOT/bin/conda" ]; then
  echo "==== installing miniconda (py37) into $CONDA_ROOT ===="
  INSTALLER=Miniconda3-py37_23.1.0-1-Linux-x86_64.sh
  TMP=$(mktemp -d)
  trap 'rm -rf "$TMP"' EXIT
  curl -fsSL -o "$TMP/$INSTALLER" \
    "https://repo.anaconda.com/miniconda/$INSTALLER"
  bash "$TMP/$INSTALLER" -b -p "$CONDA_ROOT"
fi

# ---- the environment ---------------------------------------------------------
if [ ! -d "$CONDA_ROOT/envs/$ENV_NAME" ]; then
  echo "==== creating conda env '$ENV_NAME' (python 3.7) ===="
  "$CONDA_ROOT/bin/conda" create -y -n "$ENV_NAME" python=3.7
fi

echo "==== installing tensorsignatures ===="
"$ENV_PY" -m pip install --upgrade "pip<22"
"$ENV_PY" -m pip install tensorsignatures

# Required, not optional: pip resolves protobuf to 4.x, which TensorFlow 1.15
# cannot load ("Descriptors cannot not be created directly"). 3.20.3 fixes it,
# and it must be installed AFTER tensorsignatures or its resolver undoes this.
echo "==== pinning protobuf < 3.21 (TensorFlow 1.15 cannot load 4.x) ===="
"$ENV_PY" -m pip install "protobuf<3.21"

# ---- verify ------------------------------------------------------------------
echo "==== verifying ===="
"$ENV_PY" -c "import tensorsignatures as ts, tensorflow as tf, numpy as np; \
print('tensorsignatures', ts.__version__, '| tensorflow', tf.__version__, \
'| numpy', np.__version__)"

cat <<EOF

Environment ready:  $ENV_PY

Versions used for the results reported in the paper:
  tensorsignatures 0.5.0 | tensorflow 1.15.0 | numpy 1.16.6 | protobuf 3.20.3

Next:  run_ts_sweep.sh
EOF
