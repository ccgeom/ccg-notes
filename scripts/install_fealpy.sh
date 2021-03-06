#!/usr/bin/env bash

if [ ! -f .py/lib/python${PYV}/site-packages/fealpy ]; then
    echo 'FEALPy was not installed. Installing...'
    git submodule init
    git submodule update
    cd vendors/fealpy
    pwd
    python3 setup_linux.py develop
    cd ../..
fi
