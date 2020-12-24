#!/bin/bash

# Replace $input contents with the path toward your WCSim input file
# (Note: make sure the WCSim path in RunAtStart matches the WCSim version you used to produce the file)

input=$WCSIMDIR/wcsim_hybrid.root

# test.root doesn't exist, this is just a placeholder

if [ ! -f "$input" ]; then
	echo "Set input file doesn't exist! Add your WCSim input file here. (Note: make sure the WCSim path in RunAtStart matches the WCSim version you used to produce the file)"
else
	echo "Processing ${input}..."
	./analysis $input output.root 0 0
fi

