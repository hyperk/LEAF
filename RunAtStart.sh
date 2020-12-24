#!/bin/bash

export LEAFDIR=$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LEAFDIR/lib

export WCSIMDIR=${LEAFDIR}/libWCSIM
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WCSIMDIR

############################
# ROOT/CERN
############################
# Note: official ROOT 5.34 installation on sukap is not compatible with cmake
#  Temporary use Guillaume Pronost's version. Should change in the future.
#
export ROOT_DIR=/home/pronost/software/root-5.34.38-build
# Note: There is no official ROOT 6 installation 
#  Temporary use Guillaume Pronost's version. Should change in the future.
#export ROOT_DIR=/home/pronost/software/root-6.22.00-build
echo "Load ROOT from: ${ROOT_DIR}" 
source ${ROOT_DIR}/bin/thisroot.sh

alias root='root -l'

############################
# BONSAI
############################

#Temporary, any official location?
export BONSAIDIR=/disk02/usr6/pronost/software/BonsaiHK
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BONSAIDIR

