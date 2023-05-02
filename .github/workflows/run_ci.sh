#!/usr/bin/env bash

set -e

sudo apt-get -qq install \
    gcc \
    libblas-dev \
    cmake \
    curl

python -m pip install --upgrade pip
python -m pip install pyscf
python_path=$(which python)
echo "Python path in run_ci.sh is ${python_path}"
pyscf_path=$(python -c "import pyscf; print (list(pyscf.__path__)[0], end='')")
echo "PySCF path in run_ci.sh is ${pyscf_path}"
site_packages=$(echo "${pyscf_path%/*}")
echo "site_packages is ${site_packages}"
all_python_path=$(python -c "import sys; print (sys.path)")
echo "Python paths are ${all_python_path}"
python -m pip install .

# TODO: check if pyscf code is changed using dist-info file
python -m pip uninstall -y pyscf-forge
python -m pip uninstall -y pyscf
