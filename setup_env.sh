#!/usr/bin/env bash
set -euo pipefail

#######################
# User-configurable variables
#######################
ENV_NAME="SARS-CoV-2_Deletion_detection_from_amplicon-seq"
PYTHON_VERSION="3.10"

# Channels (order matters)
CHANNEL_R="r"
CHANNEL_BIOCONDA="bioconda"
CHANNEL_CONDA_FORGE="conda-forge"

# Tool versions
STAR_VERSION="2.7.3a"
SAMTOOLS_VERSION="1.10"

#######################
# Script begins
#######################

echo ">>> Creating conda environment: ${ENV_NAME} with Python ${PYTHON_VERSION}"
conda create -y -n "${ENV_NAME}" "python=${PYTHON_VERSION}"

# Depending on your shell, you may need: source ~/miniconda3/etc/profile.d/conda.sh
echo ">>> Activating environment: ${ENV_NAME}"
conda activate "${ENV_NAME}"

echo ">>> Adding channels"
conda config --add channels "${CHANNEL_R}"
conda config --add channels "${CHANNEL_BIOCONDA}"
conda config --add channels "${CHANNEL_CONDA_FORGE}"

echo ">>> Installing mamba"
conda install -y -c conda-forge mamba

echo ">>> Installing STAR and samtools via mamba"
mamba install -y "star=${STAR_VERSION}" "samtools=${SAMTOOLS_VERSION}"

echo ">>> Installing bedtools"
conda install -y -c bioconda bedtools

echo ">>> Installing Python libraries"
conda install -y pysam pandas numpy

echo ">>> Installing Java, Perl, and Perl JSON"
conda install -y -c conda-forge openjdk perl perl-json

echo ">>> Verifying installation"
python --version
star --version
samtools --version | head -1
bedtools --version

echo ">>> Done. Environment '${ENV_NAME}' is ready."
