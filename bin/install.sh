#!/bin/bash

set -e

PREFIX=${HOME}/miniconda3

VIVI_ENV_NAME=${1-vivi}
OUTPUT=${2-/dev/stdout}

install_conda () {
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PREFIX} >> ${OUTPUT}
    export PATH=${PATH}:${PREFIX}/bin
    command -v conda > /dev/null 2>&1 || { echo "Conda still is not on the path, try installing manually"; exit 1; }
    conda update -n root conda --yes
    rm Miniconda3-latest-Linux-x86_64.sh
}

command -v conda > /dev/null 2>&1 || { echo "Conda not installed, installing now ..."; install_conda; }

conda config --prepend channels 'bushmanlab'
conda config --prepend channels 'conda-forge'
conda config --prepend channels 'r'
conda config --prepend channels 'bioconda'

# Create enviroment if it does not exist
conda env list | grep -Fxq ${VIVI_ENV_NAME} || {
    conda env create --name ${VIVI_ENV_NAME} --file bin/build.v0.1.0.yml >> ${OUTPUT}
    source activate ${VIVI_ENV_NAME}
    Rscript bin/setup.R >> ${OUTPUT}
    cd tools
    git clone https://github.com/cnobles/dualDemultiplexR.git >> ${OUTPUT}
    git clone https://github.com/cnobles/seqTrimR.git >> ${OUTPUT}
    git clone https://github.com/cnobles/seqFiltR.git >> ${OUTPUT}
    git clone https://github.com/cnobles/seqConsolidateR.git >> ${OUTPUT}
    cd ../
    echo -e "Vivi successfully installed.\n" ;
}

echo -e "To get started, ensure ${PREFIX}/bin is in your path and\n" \
  "run 'source activate ${VIVI_ENV_NAME}'\n\n" \
  "To ensure ${PREFIX}/bin is in your path each time you log in,\n" \
  "append the following to your .bashrc or .bash_profile:\n\n" \
  "# Append miniconda3/bin to path\n" \
  "export PATH='~/miniconda3/bin:${PATH}'\n"
