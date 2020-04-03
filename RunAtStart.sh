#!/bin/bash

export LEAFDIR=$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LEAFDIR/lib

export WCSIMDIR=${LEAFDIR}/libWCSIM
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WCSIMDIR

#Temporary, any official location?
export BONSAIDIR=/disk02/usr6/pronost/software/BonsaiHK
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BONSAIDIR

export ROOTSYS=/usr/local/sklib_gcc4.8.5/root_v5.34.36
source ${ROOTSYS}/bin/thisroot.sh

