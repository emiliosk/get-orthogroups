#!/usr/bin/env bash
set -e

echo "=== Treerecs Installation Helper ==="

# 1. Check conda
if ! command -v conda &> /dev/null; then
    echo "[ERROR] 'conda' was not found in PATH. Please ensure Conda or Mamba is available."
    exit 1
fi

OS="$(uname -s)"
ARCH="$(uname -m)"

echo "Detected Platform: ${OS} (${ARCH})"

if [ "$OS" = "Linux" ] || ([ "$OS" = "Darwin" ] && [ "$ARCH" = "x86_64" ]); then
    echo "Installing treerecs directly via Bioconda into active environment..."
    conda install -c bioconda treerecs -y
    echo "[OK] Treerecs installed successfully."
elif [ "$OS" = "Darwin" ] && [ "$ARCH" = "arm64" ]; then
    echo "Apple Silicon (arm64) detected."
    echo "Bioconda provides 'treerecs' compiled for osx-64 (Intel)."
    
    # Check Rosetta 2
    if ! arch -x86_64 /usr/bin/true &> /dev/null; then
        echo "[INFO] Rosetta 2 is required to run Intel binaries. Installing Rosetta..."
        softwareupdate --install-rosetta --agree-to-license
    fi

    CONDA_BASE="$(conda info --base)"
    ENV_NAME="treerecs_x86"
    TARGET_ENV="${CONDA_BASE}/envs/${ENV_NAME}"

    if [ ! -d "$TARGET_ENV" ]; then
        echo "Creating x86_64 environment '${ENV_NAME}' with osx-64 emulation..."
        CONDA_SUBDIR=osx-64 conda create -n "$ENV_NAME" -c bioconda treerecs -y
    else
        echo "Existing x86_64 environment found at ${TARGET_ENV}."
    fi

    TREERECS_EXE="${TARGET_ENV}/bin/treerecs"

    if [ ! -f "$TREERECS_EXE" ]; then
        echo "[ERROR] Treerecs binary not found at ${TREERECS_EXE}."
        exit 1
    fi

    # If inside an active conda environment, symlink into its bin/
    if [ -n "$CONDA_PREFIX" ]; then
        echo "Symlinking treerecs into active conda environment: ${CONDA_PREFIX}/bin/"
        mkdir -p "${CONDA_PREFIX}/bin"
        ln -sf "$TREERECS_EXE" "${CONDA_PREFIX}/bin/treerecs"
        echo "[OK] Symlinked to ${CONDA_PREFIX}/bin/treerecs"
    else
        echo "[NOTE] No active conda environment detected (\$CONDA_PREFIX is empty)."
        echo "To make treerecs available, run this script with your conda environment active, or set:"
        echo "  export TREERECS_BIN=\"${TREERECS_EXE}\""
    fi

    # Verify execution
    echo "Verifying treerecs execution..."
    "$TREERECS_EXE" --version
    echo "[SUCCESS] Treerecs is installed and operational on Apple Silicon via Rosetta 2!"
else
    echo "[ERROR] Unsupported platform: ${OS} ${ARCH}. Please compile treerecs from source."
    exit 1
fi
