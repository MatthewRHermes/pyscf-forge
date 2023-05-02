#!/usr/bin/env bash

set -e

sudo apt-get -qq install \
    gcc \
    libblas-dev \
    cmake \
    curl

python -m pip install --upgrade pip23.0
pip23.0 install .

# TODO: check if pyscf code is changed using dist-info file
pip23.0 uninstall -y pyscf-forge
