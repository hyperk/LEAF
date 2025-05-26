#!/bin/bash +x

# This script is used to set up the environment for running LEAF

export HitTimeLimitsNegative=-4
export HitTimeLimitsPositive=8
export SearchVtxStep=300
export SearchVtxTolerance=20
export maxHitsAngle=190
export N_Neighbors=5
export maxDistanceToNeighbors=6800

export STimePDFLimitsQueueNegative=-3;
export STimePDFLimitsQueuePositive=4;
export STimePDFLimitsQueueNegative_fullTimeWindow=0;
export STimePDFLimitsQueuePositive_fullTimeWindow=0;
export TimeWindowSizeFull=1500;
export Averaging=20; # Number of points we used to average the vertex on.

export MinimizeLimitsNegative=-700;
export MinimizeLimitsPositive=1000;

export IntegrationTimeWindow=50;

export DirTolerance=30;

export StepByStep=false; # Step by step mode was not test and is not supported by multithreading
export UseDirectionality=false;
export HighEnergy=false;
export DoubleNLL=true;
export Limit_mPMT=true;
export DirTakeAll=false;