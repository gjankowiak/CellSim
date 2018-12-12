#!/bin/sh

export CODE_BASE_DIR=/home/gaspard/Documents/git/code

export CELLSIM_DIR=${CODE_BASE_DIR}/CellSim
export VISIBILITY_DIR=${CODE_BASE_DIR}/Visibility.jl
export JULIA_TOOLS_DIR=${CODE_BASE_DIR}/julia_tools

if [ "$1" == "1.0" ];
then
    echo Processing $CELLSIM_DIR
    cd $CELLSIM_DIR && git checkout julia1.0

    echo Processing $VISIBILITY_DIR
    cd $VISIBILITY_DIR && git checkout master

    echo Processing $JULIA_TOOLS_DIR
    cd $JULIA_TOOLS_DIR && git checkout master

    exit $?
fi

if [ "$1" == "0.6" ];
then
    echo Processing $CELLSIM_DIR
    cd $CELLSIM_DIR && git checkout master

    echo Processing $VISIBILITY_DIR
    cd $VISIBILITY_DIR && git checkout julia0.6

    echo Processing $JULIA_TOOLS_DIR
    cd $JULIA_TOOLS_DIR && git checkout julia0.6
else
    echo "Usage: $0 {0.6,1.0}$"
fi
