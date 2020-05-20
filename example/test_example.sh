#!/bin/bash

cd ..
source ./RunAtStart.sh
cd example

#./analysis /disk01/usr5/bquilain/wcsim_hkhybridmpmt10pc14374100Hz_e10_center_nominal_fulltank_0hitstrigger_10000.root output.root 0 0
./analysis /disk01/usr5/bquilain/wcsim_hk20pc_4.2kHzbl_e3_center_nominal_closewall_0hitstrigger_100000.root output.root 0 0
