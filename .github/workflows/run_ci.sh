#!/usr/bin/env bash

set -e

sudo apt-get -qq install \
    gcc \
    libblas-dev \
    cmake \
    curl

python -m pip install --upgrade pip
pip install pyscf
python_path=$(which python)
echo "Python path in run_ci.sh is ${python_path}"
pyscf_path=$(python -c "import pyscf; print (list(pyscf.__path__)[0], end='')")
echo "PySCF path in run_ci.sh is ${pyscf_path}"
pip install .

# TODO: check if pyscf code is changed using dist-info file
pip uninstall -y pyscf-forge
pip uninstall -y pyscf
