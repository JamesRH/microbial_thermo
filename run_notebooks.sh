#!/usr/bin/env bash
# Launch JupyterLab on the notebooks directory, in the registered environment.
#
#   ./run_notebooks.sh              # JupyterLab
#   ./run_notebooks.sh --classic    # the classic Notebook interface
#
# Unlike setup.sh this is executed, not sourced: it activates the environment
# for its own process only and leaves your shell alone.

set -euo pipefail

env_name="microbial-thermo"
root="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export MAMBA_EXE="${MAMBA_EXE:-$HOME/miniforge3/bin/mamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/miniforge3}"

if [ ! -x "$MAMBA_EXE" ]; then
    echo "mamba not found at $MAMBA_EXE" >&2
    exit 1
fi

__mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2>/dev/null)"
eval "$__mamba_setup"
unset __mamba_setup

if ! mamba env list | grep -qE "^\s*${env_name}\s"; then
    echo "The $env_name environment does not exist yet. Run:" >&2
    echo "    source setup.sh" >&2
    exit 1
fi

mamba activate "$env_name"

# Make sure the kernel this environment provides is the one notebooks will find.
if ! jupyter kernelspec list 2>/dev/null | grep -q "$env_name"; then
    python -m ipykernel install --user \
        --name "$env_name" --display-name "$env_name" >/dev/null 2>&1
fi

interface="lab"
if [ "${1:-}" = "--classic" ]; then
    interface="notebook"
fi

if ! python -c "import jupyterlab" >/dev/null 2>&1 && [ "$interface" = "lab" ]; then
    echo "JupyterLab is not installed in $env_name; installing it..."
    mamba install -y -c conda-forge jupyterlab >/dev/null
fi

echo "Starting Jupyter ($interface) with the $env_name kernel..."
exec jupyter "$interface" --notebook-dir="$root/notebooks"
