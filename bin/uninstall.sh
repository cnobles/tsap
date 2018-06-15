#!bin/bash
# Remove created environment for vivi
VIVI_ENV_NAME=${1-vivi}

source deactivate
conda env remove -n ${VIVI_ENV_NAME} --yes
conda config --remove channels 'bushmanlab'
conda config --remove channels 'bioconda'
conda config --remove channels 'r'
conda config --remove channels 'conda-forge'
