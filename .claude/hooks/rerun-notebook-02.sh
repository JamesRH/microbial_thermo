#!/usr/bin/env bash
# Claude Code PreToolUse hook (Bash): before any `git commit`, re-run
# notebook 02 so notebooks/tower_interactive.html matches the code being
# committed, then stage the notebook pair and the page. Exit 2 blocks the
# commit and hands the error back to Claude.
set -u

command=$(jq -r '.tool_input.command // empty')
if ! grep -Eq '(^|[;&|(`[:space:]])git([[:space:]]+-C[[:space:]]+[^[:space:]]+)?[[:space:]]+commit([[:space:]]|$)' <<<"$command"; then
    exit 0
fi

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/miniforge3}"
mamba="${MAMBA_EXE:-$MAMBA_ROOT_PREFIX/bin/mamba}"
nb=notebooks/02_redox_tower_and_energy

cd "${CLAUDE_PROJECT_DIR:?}" || exit 2
log=$(mktemp)
if ! {
    "$mamba" run -n microbial-thermo jupytext --sync "$nb.ipynb" &&
        "$mamba" run -n microbial-thermo jupyter nbconvert --to notebook --execute --inplace \
            --ExecutePreprocessor.kernel_name=microbial-thermo \
            --ExecutePreprocessor.timeout=600 \
            --ExecutePreprocessor.record_timing=False "$nb.ipynb" &&
        # The saved widget state gets fresh random IDs on every run; drop it
        # so each commit does not carry a few hundred lines of churn.
        "$mamba" run -n microbial-thermo python -c '
import sys, nbformat
nb = nbformat.read(sys.argv[1], as_version=4)
nb.metadata.pop("widgets", None)
nbformat.write(nb, sys.argv[1])' "$nb.ipynb"
} >"$log" 2>&1; then
    echo "Notebook 02 failed to run, so the commit was blocked:" >&2
    tail -n 40 "$log" >&2
    rm -f "$log"
    exit 2
fi
rm -f "$log"

git add "$nb.ipynb" "$nb.py" notebooks/tower_interactive.html || exit 2
