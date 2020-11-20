#!/bin/bash

cd ..
source ./RunAtStart.sh
cd example

#./analysis /disk01/usr5/bquilain/wcsim_hkhybridmpmt10pc14374100Hz_e10_center_nominal_fulltank_0hitstrigger_10000.root output.root 0 0
#./analysis /disk01/usr5/bquilain/wcsim_hk20pc_4.2kHzbl_e3_center_nominal_closewall_0hitstrigger_100000.root output.root 0 0

# Replace ./test.root by the path toward your WCSim input file
# (Note: make sure the WCSim path in RunAtStart matches the WCSim version you used to produce the file)
input=./test.root
# test.root doesn't exist, this is just a placeholder

if [ -z $input ]; then
	echo "Set input file doesn't exist! Add your WCSim input file here. (Note: make sure the WCSim path in RunAtStart matches the WCSim version you used to produce the file)"
else
	echo "Processing ${input}..."
	./analysis $input output.root 0 0
fi

