#!/usr/bin/env bash
# Build and activate the microbial-thermo environment.
#
#   source setup.sh          # create if missing, then activate
#   source setup.sh --force  # rebuild from scratch, then activate
#
# Must be SOURCED, not executed: activating an environment changes the current
# shell, and a subshell would drop that on exit.

# --- refuse to run in a subshell --------------------------------------------
if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
    echo "setup.sh must be sourced, not executed:" >&2
    echo "    source setup.sh" >&2
    exit 1
fi

_mt_env_name="microbial-thermo"
_mt_root="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# --- put mamba on PATH ------------------------------------------------------
export MAMBA_EXE="${MAMBA_EXE:-$HOME/miniforge3/bin/mamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/miniforge3}"

if [ ! -x "$MAMBA_EXE" ]; then
    echo "mamba not found at $MAMBA_EXE" >&2
    echo "Install miniforge, or set MAMBA_EXE to your mamba binary." >&2
    return 1
fi

__mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2>/dev/null)"
if [ $? -eq 0 ]; then
    eval "$__mamba_setup"
else
    alias mamba="$MAMBA_EXE"
fi
unset __mamba_setup

# --- rebuild on request -----------------------------------------------------
if [ "$1" = "--force" ] || [ "$1" = "-f" ]; then
    echo "Removing the existing $_mt_env_name environment..."
    mamba env remove -n "$_mt_env_name" -y >/dev/null 2>&1
fi

# --- create if missing ------------------------------------------------------
if ! mamba env list | grep -qE "^\s*${_mt_env_name}\s"; then
    echo "Creating the $_mt_env_name environment (a few minutes the first time)..."
    mamba env create -f "$_mt_root/environment.yaml" || return 1
else
    echo "Environment $_mt_env_name already exists; activating it."
fi

mamba activate "$_mt_env_name" || return 1

# --- install the package itself, editable -----------------------------------
if ! python -c "import microbial_thermo" >/dev/null 2>&1; then
    echo "Installing microbial_thermo in editable mode..."
    pip install -e "$_mt_root" --quiet --no-deps || return 1
fi

# --- register the Jupyter kernel --------------------------------------------
if ! jupyter kernelspec list 2>/dev/null | grep -q "$_mt_env_name"; then
    echo "Registering the Jupyter kernel..."
    python -m ipykernel install --user \
        --name "$_mt_env_name" --display-name "$_mt_env_name" >/dev/null 2>&1
fi

echo
echo "Ready. Environment: $_mt_env_name"
python -c "import microbial_thermo as mt; print('microbial_thermo', mt.__version__, '| pygcc', mt.versions()['pygcc'])" 2>/dev/null
echo "Try:  mthermo couple HS- SO4-2"
echo "      ./run_notebooks.sh"

unset _mt_env_name _mt_root
