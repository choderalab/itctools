#! /usr/bin/env bash

if [[ ! "$TRAVIS" = true ]]; then
    export python=2.7
    export CONDA_PY=27
fi

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh  -O  Miniconda_installer.sh
elif [[ "$OSTYPE" == "darwin"* ]]; then
    wget http://repo.continuum.io/miniconda/Miniconda3-3.7.0-MacOSX-x86_64.sh -O Miniconda_installer.sh
else
    echo "Unable to detect Linux or OSX system!"
    exit 1
fi

bash Miniconda_installer.sh -b
PIP_ARGS="-U"
export PATH=$HOME/miniconda/bin:$PATH
conda update --yes conda
conda config --add channels http://conda.binstar.org/omnia
conda create --yes -n ${python} python=${python} --file requirements-conda.txt
source ${HOME}/miniconda/bin/activate ${python}

pip install -r requirements-pip.txt

python setup.py install

if [[ ! "$TRAVIS" = true ]]; then
    cd tests
    behave
fi

