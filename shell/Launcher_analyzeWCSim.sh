#!/bin/zsh +x
# Sampele shell scripts showing how to process in queue system in sukap
# This assums to be used with generateShell_analyzeWCSim.c
#
neventPerFile=1000
for e in 3 4 5 6 8 10 15;do
    for i in {0..99999};do
	f=$((${i}/$neventPerFile))
	q=$(($i - ${f}*${neventPerFile}))
	if [ $q -eq 0 ];
	then
	    echo $i
	    qsub -q all -e ${LEAFDIR}/jobs/analyzeLE_e${e}_start${i}.err -o ${LEAFDIR}/jobs/analyzeLE_e${e}_start${i}.log ${LEAFDIR}/jobs/analyzeLE_e${e}_start${i}.sh 
	fi
    done
done
