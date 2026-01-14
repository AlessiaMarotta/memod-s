#!/usr/bin/env bash
set -euo pipefail

DORADO_NAME="dorado-1.3.0-linux-x64"
DORADO_TAR="${DORADO_NAME}.tar.gz"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: No conda environment activated."
    echo "Activate a conda environment before running this script."
    exit 1
fi

INSTALL_DIR="${CONDA_PREFIX}/opt"
BIN_DIR="${INSTALL_DIR}/${DORADO_NAME}/bin"

mkdir -p "${INSTALL_DIR}"
cd "${INSTALL_DIR}"

if [[ ! -d "${DORADO_NAME}" ]]; then
    echo "Dorado installation may require several minutes due to the size of the package."
    wget https://cdn.oxfordnanoportal.com/software/analysis/${DORADO_TAR}
    tar -xzf "${DORADO_TAR}"
    rm "${DORADO_TAR}"
else
    echo "Dorado already installed at ${INSTALL_DIR}/${DORADO_NAME}"
fi

echo "${BIN_DIR}/dorado"
