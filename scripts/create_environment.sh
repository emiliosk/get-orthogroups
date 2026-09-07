#!/usr/bin/env bash
set -e

echo "=== Ensembl Orthology Pipeline Environment Setup ==="

# Check conda / mamba
CONDA_CMD=""
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
else
    echo "[ERROR] Neither 'conda' nor 'mamba' was found in PATH. Please install Conda/Mamba first."
    exit 1
fi

OS="$(uname -s)"
ARCH="$(uname -m)"

echo "Detected System: ${OS} (${ARCH})"

if [ "$OS" = "Darwin" ] && [ "$ARCH" = "arm64" ]; then
    echo "Apple Silicon detected (arm64)."
    echo "Configuring environment with osx-64 emulation via Rosetta 2 (required for treerecs)..."
    
    # Check Rosetta 2
    if ! arch -x86_64 /usr/bin/true &> /dev/null; then
        echo "[INFO] Installing Rosetta 2..."
        softwareupdate --install-rosetta --agree-to-license
    fi

    CONDA_SUBDIR=osx-64 $CONDA_CMD env create -f envs/orthology_env.yaml
    
    # Set default subdir for this environment so future conda installs also use osx-64
    CONDA_BASE="$($CONDA_CMD info --base)"
    $CONDA_CMD config --env --set subdir osx-64
else
    echo "Creating environment using standard conda configuration..."
    $CONDA_CMD env create -f envs/orthology_env.yaml
fi

echo ""
echo "[SUCCESS] Environment 'ensembl_orthology' created successfully!"
echo "Activate it using:"
echo "  conda activate ensembl_orthology"
