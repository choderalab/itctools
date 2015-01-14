#! /usr/bin/env bash

if [[ ! "$TRAVIS" = true ]]; then
    export python=2.7
    export CONDA_PY=27
fi

if [[ "${python:0:1}" == "3" ]]; then
  Miniconda="Miniconda3"
else
  Miniconda="Miniconda"
fi
	
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    wget http://repo.continuum.io/miniconda/${Miniconda}-3.7.0-Linux-x86_64.sh  -O  Miniconda_installer.sh
elif [[ "$OSTYPE" == "darwin"* ]]; then
    wget http://repo.continuum.io/miniconda/${Miniconda}-3.7.0-MacOSX-x86_64.sh -O Miniconda_installer.sh
else
    echo "Unable to detect Linux or OSX system!"
    exit 1
fi

bash Miniconda_installer.sh -b -p $HOME/miniconda
PIP_ARGS="-U"
export PATH=$HOME/miniconda/bin:$PATH
conda update --yes conda
conda install --yes --file devtools/ci/travis/requirements-conda.txt

pip install -r devtools/ci/travis/requirements-pip.txt

python setup.py install --record install.log

# On travis, tests are in travis.yml, otherwise test now
if [[ ! "$TRAVIS" = true ]]; then
    nosetests 
    behave
fi

