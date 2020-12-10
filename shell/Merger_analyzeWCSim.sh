#!/bin/zsh +x
for e in 10;do
#for e in 3 4 5 6 8 10 15;do
    hadd -f /disk01/usr5/bquilain/Analyze_hkhybridmpmt10pc14374100Hz_4.2kHzbl_${e}MeV.root /disk01/usr5/bquilain/Analyze_hkhybridmpmt10pc14374100Hz_4.2kHzbl_${e}MeV_tmp*.root
done
