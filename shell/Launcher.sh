#!/bin/zsh +x
neventPerFile=20
#for e in 3 4 5 6 8 10 15;do
for e in 10;do
    for i in {0..999};do
	f=$((${i}/$neventPerFile))
	q=$(($i - ${f}*${neventPerFile}))
	if [ $q -eq 0 ];
	then
	    echo $i
	    qsub -q all -e ${LEAFDIR}/jobs/fitLE_e${e}_start${i}.err -o ${LEAFDIR}/jobs/fitLE_e${e}_start${i}.log ${LEAFDIR}/jobs/fitLE_e${e}_start${i}.sh 
	fi
    done
done
