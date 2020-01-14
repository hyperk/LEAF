#!/bin/zsh +x
#sleep 10m
#TOTAL=$1
#COUNTER=0
#while [ $COUNTER -lt $TOTAL ]; do
#    echo Launch/RunMergerBNJ${i}.sh
#    qsub -eo -q all Launch/RunMergerBNJ${i}.sh
#    let COUNTER=COUNTER+1
#done
#neventPerFile=20
neventPerFile=10
#neventPerFile=100
#neventPerFile=100
#for e in 3 4 5 6 8 10 15;do
#for e in 4 5;do
#for e in 3 4 5;do
#for e in 6 8;do
#for e in 10 15;do
#for e in 6 8 10 15;do
#for e in 3;do
#for e in 3 4 5 6 8 10 15;do
#for e in 3 4 5;do
for e in 10;do
#    for i in {0..0};do
    for i in {0..999};do
#    for i in {0..10000};do
#    for i in {0..5000};do
	#    echo $i
	#    echo $neventPerFile
	f=$((${i}/$neventPerFile))
	q=$(($i - ${f}*${neventPerFile}))
	#    echo $f
	#    echo $q
	if [ $q -eq 0 ];
	then
	    echo $i
	    qsub -q all -e /disk01/usr5/bquilain/jobs2/fitLE_e${e}_start${i}.err -o /disk01/usr5/bquilain/jobs2/fitLE_e${e}_start${i}.log jobs/fitLE_e${e}_start${i}.sh 
#	    qsub -q all -e /disk01/usr5/bquilain/jobs2/fitLE_e${e}_start${i}.err -o /disk01/usr5/bquilain/jobs2/fitLE_e${e}_start${i}.log jobs/fitLE_e${e}_start${i}.sh 
#	    qsub -q all -e /disk01/usr5/bquilain/jobs2/inputs_e${e}_start${i}.err -o /disk01/usr5/bquilain/jobs2/inputs_e${e}_start${i}.log jobs/inputs_e${e}_start${i}.sh 
	    #-e /disk01/usr5/bquilain/jobs/ -o /disk01/usr5/bquilain/jobs/
	fi
    done
done
