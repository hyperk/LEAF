#!/bin/zsh +x
# Sampele shell scripts showing how to process in queue system in sukap
# This assums to be used with Launcher_analyzeWCSim.sh and generateShell_analyzeWCSim.c
# 
# This just merge generated root-files by them

for e in 3 4 5 6 8 10 15;do
    hadd -f /disk01/usr5/bquilain/Analyze_hkhybridmpmt10pc14374100Hz_4.2kHzbl_${e}MeV.root /disk01/usr5/bquilain/Analyze_hkhybridmpmt10pc14374100Hz_4.2kHzbl_${e}MeV_tmp*.root
done
