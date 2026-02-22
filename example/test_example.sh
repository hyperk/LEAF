#!/bin/bash
# Launch single job
#
#SBATCH --job-name=leaf                       # Job name
#SBATCH --output=/sps/t2k/lperisse/Soft/leaf/logs/analysis_%j.log                  # Standard output and error log
#SBATCH --partition=htc                              # Partition choice
#SBATCH --ntasks=1                                   # Maximum number of parallel processes
#SBATCH --mem=5G                                   # Amount of memory required
#SBATCH --time=0-01:00                              # 7 days by default on htc partition



# Replace $input contents with the path toward your WCSim input file
# (Note: make sure the WCSim path in RunAtStart matches the WCSim version you used to produce the file)

input=/sps/t2k/lperisse/Soft/wcsim/results/electron/wcsim1p12p29_electron_HK_test.root
output=/sps/t2k/lperisse/Soft/wcsim/results/electron/wcsim1p12p29_electron_HK_test.leaf.root


# Argument syntaxe to run ./analysis:
#-i  Input WCSim ROOT file
#-o  Output ROOT file
# d  Dark noise frequency of PMTs in kHz
# h  Dark noise frequency of mPMTS in kHz
# s  First event to analyze
#-e  Last event to analyze
# v VERBOSE option

if [ ! -f "$input" ]; then
	echo "Set input file doesn't exist! Add your WCSim input file here. (Note: make sure the WCSim path in RunAtStart matches the WCSim version you used to produce the file)"
else
	echo "Processing ${input}..."
	cd ${LEAFDIR}/example
	./analysis -i $input  -o $output  -s 0  -e 5  -d 4.2  -h 0.0  -v
fi



