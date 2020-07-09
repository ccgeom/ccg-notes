#!/usr/bin/env bash

if [ ! -f .py/lib/python${PYV}/site-packages/fealpy ]; then
    echo 'FEALPy was not installed. Installing...'
    git submodule init
    cd vendors/fealpy
    pwd
    python3 setup_linux.py develop
fi
