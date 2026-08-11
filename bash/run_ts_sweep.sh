#!/bin/bash
# TensorSignatures rank sweep on the chromatin-state tensor.
#
# TensorSignatures has no compressive prior, so its number of signatures has to
# be chosen by an explicit sweep and an information criterion - unlike
# SignaturePPF, where the prior selects K. That difference is itself part of the
# comparison, which is why the sweep is a first-class step here rather than a
# tuning detail.
#
# Every path is derived from the environment, so the script can be launched from
# any working directory. Each rank writes to its own sub-directory (rankNN) so
# fits never overwrite each other, and ranks that already have a complete fit are
# skipped - the sweep is restartable.
#
# Usage:   bash/run_ts_sweep.sh [rank ...]        (default: 4 5 6 7 8 9 10 11 12)
#
# Long runs:
#   nohup bash/run_ts_sweep.sh > output/02_tensorsignatures/ts_sweep.log 2>&1 &
set -u

PAPER="${SIGNATUREPPF_PAPER:-$HOME/SignaturePPF-paper}"
PY="${TENSORSIG_PYTHON:-$HOME/miniconda3/envs/tensorsig/bin/python}"
SCRIPT="$PAPER/python/run_tensorsignatures.py"
TAG="${TS_TAG:-icgc_chromatin}"
BASE="$PAPER/output/02_tensorsignatures/$TAG"

RANKS=("$@")
if [ ${#RANKS[@]} -eq 0 ]; then RANKS=(4 5 6 7 8 9 10 11 12); fi

for f in "$PY" "$SCRIPT" "$BASE/snv_counts_long.tsv.gz"; do
  if [ ! -e "$f" ]; then
    echo "MISSING: $f"
    case "$f" in
      *python)     echo "  -> run bash/setup_tensorsig_env.sh" ;;
      *.tsv.gz)    echo "  -> run Rscript R/02_tensorsignatures_prepare.R" ;;
    esac
    exit 1
  fi
done

export SIGNATUREPPF_PAPER="$PAPER"

for k in "${RANKS[@]}"; do
  D=$BASE/rank$(printf %02d "$k")
  # Skip only fits that are COMPLETE under the current script. The sentinel is
  # the newest output (ts_predicted_state_sample.tsv), so fits made before that
  # export existed are treated as stale and refitted rather than silently kept.
  if [ -f "$D/ts_predicted_state_sample.tsv" ]; then
    echo "==== rank $k already fitted ($D) - skipping ===="
    continue
  fi
  if [ -f "$D/ts_fit_summary.tsv" ]; then
    echo "==== rank $k is stale (no ts_predicted_state_sample.tsv) - refitting ===="
  fi
  echo "================ rank $k -> $D  ($(date +%H:%M:%S)) ================"
  # TensorFlow 1.15 is extremely noisy on a modern system; the filter keeps the
  # log readable without hiding anything the fit actually reports.
  "$PY" "$SCRIPT" --tag "$TAG" \
      --rank "$k" --epochs 10000 --log-step 500 --seed 1 --out "$D" \
    2>&1 | grep -vE "FutureWarning|np_resource|deprecat|Instructions|WARNING|it/s\]|XLA|StreamExecutor|tensorflow/|libcuda|cuInit|^Progress"
  echo "---- rank $k done ($(date +%H:%M:%S)) ----"
done

echo "================ SWEEP COMPLETE ================"
printf "%-6s %-14s %-8s %-12s %s\n" rank logL params AIC BIC
for d in "$BASE"/rank*/; do
  [ -f "$d/ts_fit_summary.tsv" ] || continue
  tail -n +2 "$d/ts_fit_summary.tsv" |
    awk '{printf "%-6s %-14.1f %-8s %-12.0f %.0f\n", $1, $2, $3, $5, $6}'
done

echo
echo "Next: Rscript R/03_tensorsignatures_compare.R"
