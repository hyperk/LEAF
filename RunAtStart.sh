#!/bin/bash

export LEAFDIR=$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LEAFDIR/lib

#export WCSIMDIR=${LEAFDIR}/libWCSIM
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WCSIMDIR

############################
# ROOT/CERN
############################
# Note: official ROOT 5.34 installation on sukap is not compatible with cmake
#  Temporary use Guillaume Pronost's version. Should change in the future.
#
export ROOT_DIR=/disk02/usr6/pronost/software
source ${ROOT_DIR}/root-5.34.38-build/bin/thisroot.sh

alias root='root -l'

############################
# BONSAI
############################

#Temporary, any official location?
export BONSAIDIR=/disk02/usr6/pronost/software/BonsaiHK
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BONSAIDIR

